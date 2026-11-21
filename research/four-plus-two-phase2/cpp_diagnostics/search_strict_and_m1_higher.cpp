// Search higher degrees over C via Unit-Sphere Levenberg-Marquardt:
// 1. m = 0 (strict CoreCertificate: N.grad(H) = A*D and D*div(N) - N.grad(D) = 0 — proved equivalent for ALL H \nmid D!)
//    Tests deg(D) = 7, 8, 9, 10 in pure parity!
// 2. m = 1 (D = H * D0: M.grad(H) = 0 and 3*(D0*div(M) - M.grad(D0)) + 2*A*D0^2 = 0)
//    Tests deg(D0) = 3, 4, 5, 6 in pure parity (i.e. deg(D) = 9, 10, 11, 12)!

#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <vector>

using cd = std::complex<double>;

struct Mon4 {
    int e[4];
    bool operator<(const Mon4& o) const {
        for (int i=0;i<4;++i) if (e[i]!=o.e[i]) return e[i]<o.e[i];
        return false;
    }
    int pa() const { return (e[0] + e[1]) & 1; }
    int pc() const { return (e[2] + e[3]) & 1; }
};

static std::vector<Mon4> gen_mons4(int deg, int req_pa = -1, int req_pc = -1) {
    std::vector<Mon4> res;
    for (int a = deg; a >= 0; --a)
        for (int b = deg - a; b >= 0; --b)
            for (int c = deg - a - b; c >= 0; --c) {
                int d = deg - a - b - c;
                Mon4 m{{a, b, c, d}};
                if (req_pa != -1 && m.pa() != req_pa) continue;
                if (req_pc != -1 && m.pc() != req_pc) continue;
                res.push_back(m);
            }
    return res;
}

using CPoly4 = std::map<Mon4, double>;
static void cpoly_add(CPoly4& p, const Mon4& m, double c) {
    if (std::abs(c) < 1e-15) return;
    p[m] += c;
    if (std::abs(p[m]) < 1e-15) p.erase(m);
}
static CPoly4 cpoly_mul(const CPoly4& A, const CPoly4& B) {
    CPoly4 C;
    for (auto& [ma, ca] : A)
        for (auto& [mb, cb] : B)
            cpoly_add(C, Mon4{{ma.e[0]+mb.e[0], ma.e[1]+mb.e[1], ma.e[2]+mb.e[2], ma.e[3]+mb.e[3]}}, ca * cb);
    return C;
}
static CPoly4 cpoly_deriv(const CPoly4& A, int v) {
    CPoly4 C;
    for (auto& [m, c] : A) if (m.e[v] > 0) {
        Mon4 md = m; double ev = md.e[v]; md.e[v]--;
        cpoly_add(C, md, c * ev);
    }
    return C;
}

static std::vector<cd> solve_hermitian(std::vector<std::vector<cd>> H, std::vector<cd> b) {
    int n = (int)b.size();
    for (int i = 0; i < n; ++i) {
        int piv = i; double best = std::abs(H[i][i]);
        for (int r = i + 1; r < n; ++r) { double v = std::abs(H[r][i]); if (v > best) { best = v; piv = r; } }
        if (best < 1e-22) continue;
        if (piv != i) { std::swap(H[i], H[piv]); std::swap(b[i], b[piv]); }
        cd inv = cd(1.0, 0.0) / H[i][i];
        for (int j = i; j < n; ++j) H[i][j] *= inv;
        b[i] *= inv;
        for (int r = 0; r < n; ++r) if (r != i && std::abs(H[r][i]) > 0.0) {
            cd f = H[r][i];
            for (int j = i; j < n; ++j) H[r][j] -= f * H[i][j];
            b[r] -= f * b[i];
        }
    }
    return b;
}

struct BilinearTerm { int r, i, j; double val; };

static void run_lm_bilinear(const std::string& label, int nD, int nK, int m_eq, const std::vector<BilinearTerm>& terms, int max_it = 180) {
    int nVars = nD + nK;
    std::mt19937_64 rng(777);
    std::normal_distribution<double> nd(0.0, 1.0);

    double best_min = 1e99;
    for (int st = 0; st < 2; ++st) {
        std::vector<cd> x(nVars), d0_ref(nD);
        double nrmD = 0.0;
        for (int i=0;i<nD;++i) { x[i] = cd(nd(rng), nd(rng)); nrmD += std::norm(x[i]); }
        nrmD = std::sqrt(nrmD);
        for (int i=0;i<nD;++i) { x[i] /= nrmD; d0_ref[i] = x[i]; }
        for (int i=nD;i<nVars;++i) x[i] = cd(0,0);

        double lam = 1e-3;
        auto eval_poly = [&](const std::vector<cd>& xv, std::vector<cd>& res) {
            res.assign(m_eq + 1, cd(0,0));
            for (const auto& t : terms) res[t.r] += t.val * xv[t.i] * xv[t.j];
            cd d(0,0);
            for (int i=0;i<nD;++i) d += std::conj(d0_ref[i]) * xv[i];
            res[m_eq] = (d - cd(1.0, 0.0)) * 10.0;
            double nr2 = 0.0;
            for (int r=0;r<=m_eq;++r) nr2 += std::norm(res[r]);
            return std::sqrt(nr2);
        };

        std::vector<cd> res;
        double cur_err = eval_poly(x, res);
        for (int it = 1; it <= max_it; ++it) {
            std::vector<std::vector<cd>> J(m_eq + 1, std::vector<cd>(nVars, cd(0,0)));
            for (const auto& t : terms) {
                if (t.i == t.j) J[t.r][t.i] += (2.0 * t.val) * x[t.i];
                else { J[t.r][t.i] += t.val * x[t.j]; J[t.r][t.j] += t.val * x[t.i]; }
            }
            for (int i=0;i<nD;++i) J[m_eq][i] = std::conj(d0_ref[i]) * 10.0;

            std::vector<std::vector<cd>> H(nVars, std::vector<cd>(nVars, cd(0,0)));
            std::vector<cd> g(nVars, cd(0,0));
            for (int rr=0;rr<=m_eq;++rr) {
                for (int i=0;i<nVars;++i) {
                    if (J[rr][i] == cd(0,0)) continue;
                    cd jc = std::conj(J[rr][i]);
                    g[i] -= jc * res[rr];
                    for (int j=i;j<nVars;++j) H[i][j] += jc * J[rr][j];
                }
            }
            for (int i=0;i<nVars;++i) {
                for (int j=0;j<i;++j) H[i][j] = std::conj(H[j][i]);
                H[i][i] += lam;
            }
            auto dx = solve_hermitian(std::move(H), std::move(g));
            std::vector<cd> xn(nVars);
            for (int i=0;i<nVars;++i) xn[i] = x[i] + dx[i];
            std::vector<cd> rn;
            double nerr = eval_poly(xn, rn);
            if (nerr < cur_err) {
                x.swap(xn); res.swap(rn); cur_err = nerr;
                lam = std::max(1e-15, lam * 0.3);
            } else {
                lam = std::min(1e8, lam * 5.0);
            }
            if (it % 45 == 0) {
                double nd2 = 0;
                for (int i=0;i<nD;++i) nd2 += std::norm(x[i]);
                nd2 = std::sqrt(nd2);
                for (int i=0;i<nVars;++i) x[i] /= nd2;
                for (int i=0;i<nD;++i) d0_ref[i] = x[i];
                cur_err = eval_poly(x, res);
            }
            if (cur_err < 1e-11) break;
        }
        double p_err2 = 0;
        for (int rr=0;rr<m_eq;++rr) p_err2 += std::norm(res[rr]);
        best_min = std::min(best_min, std::sqrt(p_err2));
    }
    std::cout << label << " (nD=" << nD << ", nK=" << nK << ", eqs=" << m_eq
              << ") -> min ||eq||_2(at ||D||=1) = " << std::scientific << std::setprecision(4) << best_min
              << (best_min < 1e-8 ? "  *** COMPLEX SOLUTION DETECTED!! ***" : "") << "\n";
}

static void test_strict_deg(int degD, int par_ab, int par_cd, const CPoly4& H7, const CPoly4& A7) {
    CPoly4 dH[4];
    for (int i=0;i<4;++i) dH[i] = cpoly_deriv(H7, i);
    auto D_mons = gen_mons4(degD, par_ab, par_cd);
    std::vector<Mon4> N_mons[4];
    for (int i=0;i<4;++i) {
        if (i < 2) N_mons[i] = gen_mons4(degD + 1, par_ab, par_cd ^ 1);
        else       N_mons[i] = gen_mons4(degD + 1, par_ab ^ 1, par_cd);
    }
    int nD_raw = (int)D_mons.size(), off[5] = {0};
    for (int i=0;i<4;++i) off[i+1] = off[i] + (int)N_mons[i].size();
    int nN_raw = off[4], tot1 = nN_raw + nD_raw;

    auto eq1 = gen_mons4(degD + 6, par_ab ^ 1, par_cd ^ 1);
    std::map<Mon4, int> eq1_idx;
    for (int i=0;i<(int)eq1.size();++i) eq1_idx[eq1[i]] = i;

    std::vector<std::vector<double>> M1(eq1.size(), std::vector<double>(tot1, 0.0));
    for (int i=0;i<4;++i)
        for (int j=0;j<(int)N_mons[i].size();++j)
            for (auto& [m, c] : cpoly_mul(CPoly4{{N_mons[i][j], 1.0}}, dH[i]))
                M1[eq1_idx.at(m)][off[i]+j] += c;
    for (int j=0;j<nD_raw;++j)
        for (auto& [m, c] : cpoly_mul(CPoly4{{D_mons[j], 1.0}}, A7))
            M1[eq1_idx.at(m)][nN_raw+j] -= c;

    int nrows = (int)M1.size(), r = 0;
    std::vector<int> piv_col, col_to_row(tot1, -1);
    for (int c=0;c<tot1 && r<nrows;++c) {
        int best=-1; double bv=1e-9;
        for (int i=r;i<nrows;++i) if (std::abs(M1[i][c]) > bv) { bv=std::abs(M1[i][c]); best=i; }
        if (best==-1) continue;
        std::swap(M1[r], M1[best]);
        double inv = 1.0 / M1[r][c];
        for (int j=c;j<tot1;++j) M1[r][j] *= inv;
        for (int i=0;i<nrows;++i) if (i!=r && std::abs(M1[i][c])>1e-15) {
            double f = M1[i][c];
            for (int j=c;j<tot1;++j) M1[i][j] -= f * M1[r][j];
        }
        piv_col.push_back(c); col_to_row[c] = r++;
    }
    std::vector<std::vector<double>> K_vecs, E1_vecs;
    for (int c=0;c<tot1;++c) if (col_to_row[c]==-1) {
        std::vector<double> v(tot1, 0.0); v[c]=1.0;
        for (int pr=0;pr<r;++pr) v[piv_col[pr]] = -M1[pr][c];
        bool dz=true;
        for (int j=0;j<nD_raw;++j) if (std::abs(v[nN_raw+j])>1e-9) { dz=false; break; }
        if (dz) K_vecs.push_back(std::move(v));
        else    E1_vecs.push_back(std::move(v));
    }
    int nD = (int)E1_vecs.size(), nK = (int)K_vecs.size();

    std::vector<CPoly4> D_b(nD), dD_b[4], divN0(nD), N0_b[4], K_b[4], divK(nK);
    for (int i=0;i<4;++i) { dD_b[i].resize(nD); N0_b[i].resize(nD); K_b[i].resize(nK); }
    for (int a=0;a<nD;++a) {
        for (int j=0;j<nD_raw;++j) if (std::abs(E1_vecs[a][nN_raw+j])>1e-12) D_b[a][D_mons[j]] = E1_vecs[a][nN_raw+j];
        for (int i=0;i<4;++i) {
            dD_b[i][a] = cpoly_deriv(D_b[a], i);
            for (int j=0;j<(int)N_mons[i].size();++j) if (std::abs(E1_vecs[a][off[i]+j])>1e-12) N0_b[i][a][N_mons[i][j]] = E1_vecs[a][off[i]+j];
            for (auto& [m, c] : cpoly_deriv(N0_b[i][a], i)) cpoly_add(divN0[a], m, c);
        }
    }
    for (int b=0;b<nK;++b) {
        for (int i=0;i<4;++i) {
            for (int j=0;j<(int)N_mons[i].size();++j) if (std::abs(K_vecs[b][off[i]+j])>1e-12) K_b[i][b][N_mons[i][j]] = K_vecs[b][off[i]+j];
            for (auto& [m, c] : cpoly_deriv(K_b[i][b], i)) cpoly_add(divK[b], m, c);
        }
    }

    auto div_eq = gen_mons4(2 * degD, 1, 1);
    std::map<Mon4, int> div_idx;
    for (int i=0;i<(int)div_eq.size();++i) div_idx[div_eq[i]] = i;

    std::vector<BilinearTerm> terms;
    for (int a1=0;a1<nD;++a1) for (int a2=a1;a2<nD;++a2) {
        CPoly4 P = cpoly_mul(D_b[a1], divN0[a2]);
        for (int i=0;i<4;++i) for (auto& [m,c] : cpoly_mul(N0_b[i][a2], dD_b[i][a1])) cpoly_add(P, m, -c);
        if (a1 < a2) {
            for (auto& [m,c] : cpoly_mul(D_b[a2], divN0[a1])) cpoly_add(P, m, c);
            for (int i=0;i<4;++i) for (auto& [m,c] : cpoly_mul(N0_b[i][a1], dD_b[i][a2])) cpoly_add(P, m, -c);
        }
        for (auto& [m, c] : P) terms.push_back({div_idx.at(m), a1, a2, c});
    }
    for (int a=0;a<nD;++a) for (int b=0;b<nK;++b) {
        CPoly4 P = cpoly_mul(D_b[a], divK[b]);
        for (int i=0;i<4;++i) for (auto& [m,c] : cpoly_mul(K_b[i][b], dD_b[i][a])) cpoly_add(P, m, -c);
        for (auto& [m, c] : P) terms.push_back({div_idx.at(m), a, nD + b, c});
    }
    std::string lbl = "[m=0 strict deg=" + std::to_string(degD) + " par=(" + std::to_string(par_ab) + "," + std::to_string(par_cd) + ")]";
    run_lm_bilinear(lbl, nD, nK, (int)div_eq.size(), terms, 150);
}

static void test_m1_deg(int d0, int par_ab, int par_cd, const CPoly4& H7, const CPoly4& A7) {
    CPoly4 dH[4];
    for (int i=0;i<4;++i) dH[i] = cpoly_deriv(H7, i);
    auto D0_mons = gen_mons4(d0, par_ab, par_cd);
    int nD = (int)D0_mons.size();
    std::vector<CPoly4> D0_b(nD), dD0_b[4];
    for (int i=0;i<4;++i) dD0_b[i].resize(nD);
    for (int a=0;a<nD;++a) {
        D0_b[a][D0_mons[a]] = 1.0;
        for (int i=0;i<4;++i) dD0_b[i][a] = cpoly_deriv(D0_b[a], i);
    }

    // Build independent basis of Koszul syzygies M of deg d0 + 7 (mod Kosz_2 syzygies!)
    std::vector<Mon4> M_mons[4];
    for (int i=0;i<4;++i) {
        if (i < 2) M_mons[i] = gen_mons4(d0 + 7, par_ab, par_cd ^ 1);
        else       M_mons[i] = gen_mons4(d0 + 7, par_ab ^ 1, par_cd);
    }
    int off[5] = {0};
    for (int i=0;i<4;++i) off[i+1] = off[i] + (int)M_mons[i].size();
    int nM_raw = off[4];
    auto eq1 = gen_mons4(d0 + 12, par_ab ^ 1, par_cd ^ 1);
    std::map<Mon4, int> eq1_idx;
    for (int i=0;i<(int)eq1.size();++i) eq1_idx[eq1[i]] = i;

    std::vector<std::vector<double>> M1(eq1.size(), std::vector<double>(nM_raw, 0.0));
    for (int i=0;i<4;++i)
        for (int j=0;j<(int)M_mons[i].size();++j)
            for (auto& [m, c] : cpoly_mul(CPoly4{{M_mons[i][j], 1.0}}, dH[i]))
                M1[eq1_idx.at(m)][off[i]+j] += c;

    int nrows = (int)M1.size(), r = 0;
    std::vector<int> piv_col, col_to_row(nM_raw, -1);
    for (int c=0;c<nM_raw && r<nrows;++c) {
        int best=-1; double bv=1e-9;
        for (int i=r;i<nrows;++i) if (std::abs(M1[i][c]) > bv) { bv=std::abs(M1[i][c]); best=i; }
        if (best==-1) continue;
        std::swap(M1[r], M1[best]);
        double inv = 1.0 / M1[r][c];
        for (int j=c;j<nM_raw;++j) M1[r][j] *= inv;
        for (int i=0;i<nrows;++i) if (i!=r && std::abs(M1[i][c])>1e-15) {
            double f = M1[i][c];
            for (int j=c;j<nM_raw;++j) M1[i][j] -= f * M1[r][j];
        }
        piv_col.push_back(c); col_to_row[c] = r++;
    }
    std::vector<std::vector<double>> K_vecs;
    for (int c=0;c<nM_raw;++c) if (col_to_row[c]==-1) {
        std::vector<double> v(nM_raw, 0.0); v[c]=1.0;
        for (int pr=0;pr<r;++pr) v[piv_col[pr]] = -M1[pr][c];
        K_vecs.push_back(std::move(v));
    }
    int nK = (int)K_vecs.size();

    std::vector<CPoly4> K_b[4], divK(nK);
    for (int i=0;i<4;++i) K_b[i].resize(nK);
    for (int b=0;b<nK;++b) {
        for (int i=0;i<4;++i) {
            for (int j=0;j<(int)M_mons[i].size();++j) if (std::abs(K_vecs[b][off[i]+j])>1e-12) K_b[i][b][M_mons[i][j]] = K_vecs[b][off[i]+j];
            for (auto& [m, c] : cpoly_deriv(K_b[i][b], i)) cpoly_add(divK[b], m, c);
        }
    }

    auto eq2 = gen_mons4(2 * d0 + 6, 1, 1);
    std::map<Mon4, int> eq2_idx;
    for (int i=0;i<(int)eq2.size();++i) eq2_idx[eq2[i]] = i;

    std::vector<BilinearTerm> terms;
    for (int a1=0;a1<nD;++a1) for (int a2=a1;a2<nD;++a2) {
        double mult = (a1 == a2) ? 2.0 : 4.0;
        for (auto& [m, c] : cpoly_mul(A7, cpoly_mul(D0_b[a1], D0_b[a2])))
            terms.push_back({eq2_idx.at(m), a1, a2, mult * c});
    }
    for (int a=0;a<nD;++a) for (int b=0;b<nK;++b) {
        CPoly4 P = cpoly_mul(D0_b[a], divK[b]);
        for (int i=0;i<4;++i) for (auto& [m,c] : cpoly_mul(K_b[i][b], dD0_b[i][a])) cpoly_add(P, m, -c);
        for (auto& [m, c] : P) terms.push_back({eq2_idx.at(m), a, nD + b, 3.0 * c});
    }
    std::string lbl = "[m=1 d0=" + std::to_string(d0) + " (degD=" + std::to_string(d0+6) + ") par=(" + std::to_string(par_ab) + "," + std::to_string(par_cd) + ")]";
    run_lm_bilinear(lbl, nD, nK, (int)eq2.size(), terms, 150);
}

int main() {
    CPoly4 H7, A7;
    H7[Mon4{{6,0,0,0}}] = 7; H7[Mon4{{0,6,0,0}}] = 7;
    H7[Mon4{{0,0,6,0}}] = 7; H7[Mon4{{0,0,0,6}}] = 7;
    H7[Mon4{{5,1,0,0}}] = 5; H7[Mon4{{1,5,0,0}}] = 10;
    H7[Mon4{{0,0,5,1}}] = 5; H7[Mon4{{0,0,1,5}}] = 15;
    H7[Mon4{{3,1,1,1}}] = 5; H7[Mon4{{1,1,3,1}}] = 10;
    A7[Mon4{{3,0,3,0}}] = 7;

    // m = 0 strict:
    test_strict_deg(7, 1, 0, H7, A7);
    test_strict_deg(7, 0, 1, H7, A7);
    test_strict_deg(8, 0, 0, H7, A7);
    test_strict_deg(8, 1, 1, H7, A7);
    test_strict_deg(9, 1, 0, H7, A7);
    test_strict_deg(10, 0, 0, H7, A7);

    // m = 1 (D = H * D0):
    test_m1_deg(3, 1, 0, H7, A7);
    test_m1_deg(4, 0, 0, H7, A7);
    test_m1_deg(4, 1, 1, H7, A7);
    test_m1_deg(5, 1, 0, H7, A7);
    return 0;
}
