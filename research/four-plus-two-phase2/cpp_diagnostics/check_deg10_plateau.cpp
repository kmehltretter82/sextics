// Deep convergence test on m=0 strict deg=10, 11, 12 (up to 600 iterations)
// using FAST QR/Cholesky Variable-Projection or Sparse Jacobian!
// Wait: forfixed D, N.grad(H)=A*D and D*div(N) - N.grad(D) = 0 is Q(D,D) + B(D,K) = 0!
// For ANY fixed D with ||D||_2 = 1, the optimal K minimizing ||Q(D,D) + B(D,K)||_2 is simply
// the linear least-squares residual of Q(D,D) modulo Im(B(D, .))!!
// Wait: computing min_K ||B(D) K + Q(D,D)||_2 at a single random D vs optimizing D on the unit sphere!

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

static void check_deg10_plateau(int degD, int par_ab, int par_cd, const CPoly4& H7, const CPoly4& A7) {
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
    // ORTHONORMALIZE E1_vecs on their D-coordinates so ||x_D||_2 IS LITERALLY ||D||_{L2(monomials)}!!
    // AHA!! Look! Earlier ||x_D||_2 was in the non-orthonormal RREF basis where big pivot entries inflated ||D||!!
    std::vector<CPoly4> D_b(nD), dD_b[4], divN0(nD), N0_b[4], K_b[4], divK(nK);
    for (int i=0;i<4;++i) { dD_b[i].resize(nD); N0_b[i].resize(nD); K_b[i].resize(nK); }
    for (int a=0;a<nD;++a) {
        double nrm = 0.0;
        for (int j=0;j<nD_raw;++j) nrm += E1_vecs[a][nN_raw+j] * E1_vecs[a][nN_raw+j];
        nrm = std::sqrt(nrm);
        for (double& z : E1_vecs[a]) z /= nrm;
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
    int m_eq = (int)div_eq.size();

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

    int nVars = nD + nK;
    std::mt19937_64 rng(999);
    std::normal_distribution<double> nd(0.0, 1.0);
    std::vector<cd> x(nVars), d0_ref(nD);
    double nrmD = 0.0;
    for (int i=0;i<nD;++i) { x[i] = cd(nd(rng), nd(rng)); nrmD += std::norm(x[i]); }
    nrmD = std::sqrt(nrmD);
    for (int i=0;i<nD;++i) { x[i] /= nrmD; d0_ref[i] = x[i]; }
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
    std::cout << "\n=== Normalized D-L2 check deg=" << degD << " (nD=" << nD << ", nK=" << nK << ", eqs=" << m_eq << ") ===\n";
    for (int it = 1; it <= 300; ++it) {
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
        if (nerr < cur_err) { x.swap(xn); res.swap(rn); cur_err = nerr; lam = std::max(1e-15, lam * 0.3); }
        else lam = std::min(1e8, lam * 5.0);
        if (it % 60 == 0) {
            // Compute TRUE monomial L2 norm of polynomial D = sum x[a] D_b[a]!
            std::map<Mon4, cd> Dpoly;
            for (int a=0;a<nD;++a) for (auto& [m, c] : D_b[a]) Dpoly[m] += x[a] * c;
            double true_nrmD = 0;
            for (auto& [m, z] : Dpoly) true_nrmD += std::norm(z);
            true_nrmD = std::sqrt(true_nrmD);
            for (int i=0;i<nVars;++i) x[i] /= true_nrmD;
            double xdn = 0;
            for (int i=0;i<nD;++i) xdn += std::norm(x[i]);
            for (int i=0;i<nD;++i) d0_ref[i] = x[i] / xdn;
            cur_err = eval_poly(x, res);
            double pe2 = 0;
            for (int rr=0;rr<m_eq;++rr) pe2 += std::norm(res[rr]);
            std::cout << "  it " << it << ": ||eq||_2 (at true ||D||_2=1) = " << std::scientific << std::setprecision(4) << std::sqrt(pe2) << "\n";
        }
    }
}

int main() {
    CPoly4 H7, A7;
    H7[Mon4{{6,0,0,0}}] = 7; H7[Mon4{{0,6,0,0}}] = 7;
    H7[Mon4{{0,0,6,0}}] = 7; H7[Mon4{{0,0,0,6}}] = 7;
    H7[Mon4{{5,1,0,0}}] = 5; H7[Mon4{{1,5,0,0}}] = 10;
    H7[Mon4{{0,0,5,1}}] = 5; H7[Mon4{{0,0,1,5}}] = 15;
    H7[Mon4{{3,1,1,1}}] = 5; H7[Mon4{{1,1,3,1}}] = 10;
    A7[Mon4{{3,0,3,0}}] = 7;

    check_deg10_plateau(10, 0, 0, H7, A7);
    return 0;
}
