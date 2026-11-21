// Fast exact modular source-image auditor for four-plus-two-residual:
// Computes dim(W_D) = 40 - rank([M_S M_N]) + rank(M_N) and dim(W_D / E) = dim(W_D) - 8
// across:
//   1. All 19 rational covariants & structured denominators (checked at TWO large primes p1, p2).
//   2. Positive control: uncoupled binary-sum core H_0 = f_a + f_c at p = 23 with binary-contact
//      divisors, verifying our auditor detects genuine transcendental classes in W_D / E_0!
//   3. Coupled core H at p = 23 (and F_{23^2} where needed) on benediamond's nodal-contact surfaces
//      and elliptic cones: a^3 + B(c,d), b^3 + B(c,d), H - C(c,d)^3, B_a^2 - B_c^2, etc.!

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <vector>

using u64 = uint64_t;

struct ModCtx {
    u64 P;
    u64 add(u64 a, u64 b) const { u64 s = a + b; return s >= P ? s - P : s; }
    u64 sub(u64 a, u64 b) const { return a >= b ? a - b : a + P - b; }
    u64 mul(u64 a, u64 b) const { return (unsigned __int128)a * b % P; }
    u64 pow_m(u64 a, u64 e) const {
        u64 r = 1;
        while (e) { if (e & 1) r = mul(r, a); a = mul(a, a); e >>= 1; }
        return r;
    }
    u64 inv(u64 a) const { return pow_m(a % P, P - 2); }
};

struct Mon4 {
    int e[4];
    bool operator==(const Mon4& o) const {
        return e[0]==o.e[0] && e[1]==o.e[1] && e[2]==o.e[2] && e[3]==o.e[3];
    }
    bool operator<(const Mon4& o) const {
        for (int i=0;i<4;++i) if (e[i]!=o.e[i]) return e[i]<o.e[i];
        return false;
    }
    int pa() const { return (e[0] + e[1]) & 1; }
    int pc() const { return (e[2] + e[3]) & 1; }
    int deg() const { return e[0] + e[1] + e[2] + e[3]; }
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

using Poly4 = std::map<Mon4, u64>;

static void poly_add_term(Poly4& p, const Mon4& m, u64 c, const ModCtx& M) {
    c %= M.P;
    if (!c) return;
    auto it = p.find(m);
    if (it == p.end()) p[m] = c;
    else {
        u64 s = M.add(it->second, c);
        if (s == 0) p.erase(it);
        else it->second = s;
    }
}

static Poly4 poly_add(const Poly4& A, const Poly4& B, const ModCtx& M) {
    Poly4 C = A;
    for (auto& [m, c] : B) poly_add_term(C, m, c, M);
    return C;
}

static Poly4 poly_sub(const Poly4& A, const Poly4& B, const ModCtx& M) {
    Poly4 C = A;
    for (auto& [m, c] : B) poly_add_term(C, m, M.sub(0, c), M);
    return C;
}

static Poly4 poly_scale(const Poly4& A, u64 s, const ModCtx& M) {
    Poly4 C;
    s %= M.P;
    if (!s) return C;
    for (auto& [m, c] : A) poly_add_term(C, m, M.mul(c, s), M);
    return C;
}

static Poly4 poly_mul(const Poly4& A, const Poly4& B, const ModCtx& M) {
    Poly4 C;
    for (auto& [ma, ca] : A)
        for (auto& [mb, cb] : B) {
            Mon4 mc{{ma.e[0]+mb.e[0], ma.e[1]+mb.e[1], ma.e[2]+mb.e[2], ma.e[3]+mb.e[3]}};
            poly_add_term(C, mc, M.mul(ca, cb), M);
        }
    return C;
}

static Poly4 poly_pow(Poly4 A, int e, const ModCtx& M) {
    Poly4 R{{Mon4{{0,0,0,0}}, 1}};
    while (e > 0) {
        if (e & 1) R = poly_mul(R, A, M);
        e >>= 1;
        if (e) A = poly_mul(A, A, M);
    }
    return R;
}

static Poly4 poly_deriv(const Poly4& A, int var, const ModCtx& M) {
    Poly4 C;
    for (auto& [m, c] : A) {
        if (m.e[var] > 0) {
            Mon4 md = m;
            u64 coeff = M.mul(c, (u64)md.e[var] % M.P);
            md.e[var]--;
            poly_add_term(C, md, coeff, M);
        }
    }
    return C;
}

using SparseCol = std::vector<std::pair<int, u64>>;

static void norm_col(SparseCol& col, const ModCtx& M) {
    if (col.empty()) return;
    std::sort(col.begin(), col.end(), [](const auto& a, const auto& b){ return a.first < b.first; });
    size_t w = 0;
    for (size_t i = 0; i < col.size(); ++i) {
        u64 val = col[i].second % M.P;
        if (!val) continue;
        if (w > 0 && col[w-1].first == col[i].first) {
            col[w-1].second = M.add(col[w-1].second, val);
            if (!col[w-1].second) --w;
        } else {
            col[w++] = {col[i].first, val};
        }
    }
    col.resize(w);
}

// Column-wise sparse rank: insert columns from N first, recording rank_N, then insert S columns, recording rank_full!
// Notice how inserting N columns FIRST and THEN the 40 S columns computes BOTH rank(M_N) AND rank([M_N M_S]) = rank([M_S M_N])
// IN ONE SINGLE ELIMINATION PASS!!
struct SourceAuditResult {
    int degD;
    int nRows;
    int nN;
    int rank_N;
    int rank_full;
    int dim_WD;
    int new_classes;
    bool hits_A;
};

static SourceAuditResult audit_WD(const Poly4& D, const Poly4& H, const Poly4& A_src, const ModCtx& M, int base_E = 8) {
    int degD = D.begin()->first.deg();
    int pa = D.begin()->first.pa();
    int pc = D.begin()->first.pc();
    bool pure = true;
    for (auto& [m, _] : D) {
        if (m.pa() != pa || m.pc() != pc) { pure = false; break; }
    }

    Poly4 dH[4], dD[4];
    for (int i = 0; i < 4; ++i) {
        dH[i] = poly_deriv(H, i, M);
        dD[i] = poly_deriv(D, i, M);
    }

    std::vector<Mon4> N_mons[4];
    for (int i = 0; i < 4; ++i) {
        if (!pure) N_mons[i] = gen_mons4(degD + 1, -1, -1);
        else if (i < 2) N_mons[i] = gen_mons4(degD + 1, pa, pc ^ 1);
        else            N_mons[i] = gen_mons4(degD + 1, pa ^ 1, pc);
    }
    int nN = 0;
    for (int i = 0; i < 4; ++i) nN += (int)N_mons[i].size();

    auto eq_mons = (!pure) ? gen_mons4(2 * degD + 6, -1, -1) : gen_mons4(2 * degD + 6, 1, 1);
    std::map<Mon4, int> eq_idx;
    for (int i = 0; i < (int)eq_mons.size(); ++i) eq_idx[eq_mons[i]] = i;
    int nRows = (int)eq_mons.size();

    Poly4 HD3 = poly_scale(poly_mul(H, D, M), 3, M);
    Poly4 D2 = poly_mul(D, D, M);
    Poly4 cNi[4];
    for (int i = 0; i < 4; ++i) {
        Poly4 t1 = poly_scale(poly_mul(D, dH[i], M), 2, M);
        Poly4 t2 = poly_scale(poly_mul(H, dD[i], M), 3, M);
        cNi[i] = poly_sub(Poly4{}, poly_add(t1, t2, M), M);
    }

    std::vector<SparseCol> n_cols;
    n_cols.reserve(nN);
    for (int i = 0; i < 4; ++i) {
        for (const auto& mn : N_mons[i]) {
            SparseCol col;
            for (auto& [mc, cc] : cNi[i]) {
                Mon4 mt{{mc.e[0]+mn.e[0], mc.e[1]+mn.e[1], mc.e[2]+mn.e[2], mc.e[3]+mn.e[3]}};
                col.emplace_back(eq_idx.at(mt), cc);
            }
            if (mn.e[i] > 0) {
                Mon4 dm = mn;
                u64 ev = (u64)dm.e[i] % M.P;
                dm.e[i]--;
                for (auto& [mc, cc] : HD3) {
                    Mon4 mt{{mc.e[0]+dm.e[0], mc.e[1]+dm.e[1], mc.e[2]+dm.e[2], mc.e[3]+dm.e[3]}};
                    col.emplace_back(eq_idx.at(mt), M.mul(cc, ev));
                }
            }
            norm_col(col, M);
            n_cols.push_back(std::move(col));
        }
    }

    std::vector<int> order(nN);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b){ return n_cols[a].size() < n_cols[b].size(); });

    std::vector<SparseCol> pivots(nRows);
    std::vector<bool> has_piv(nRows, false);
    SparseCol tmp;

    auto reduce_col = [&](SparseCol& cur) {
        while (!cur.empty()) {
            int lead = cur[0].first;
            if (!has_piv[lead]) return lead;
            u64 fac = cur[0].second;
            const SparseCol& piv = pivots[lead];
            tmp.clear();
            tmp.reserve(cur.size() + piv.size());
            size_t i = 1, j = 1;
            while (i < cur.size() && j < piv.size()) {
                if (cur[i].first < piv[j].first) tmp.push_back(cur[i++]);
                else if (piv[j].first < cur[i].first) {
                    u64 v = M.sub(0, M.mul(fac, piv[j].second));
                    if (v) tmp.emplace_back(piv[j].first, v);
                    ++j;
                } else {
                    u64 v = M.sub(cur[i].second, M.mul(fac, piv[j].second));
                    if (v) tmp.emplace_back(cur[i].first, v);
                    ++i; ++j;
                }
            }
            while (i < cur.size()) tmp.push_back(cur[i++]);
            while (j < piv.size()) {
                u64 v = M.sub(0, M.mul(fac, piv[j].second));
                if (v) tmp.emplace_back(piv[j].first, v);
                ++j;
            }
            cur.swap(tmp);
        }
        return -1;
    };

    int rank_N = 0;
    for (int idx : order) {
        SparseCol cur = std::move(n_cols[idx]);
        int lead = reduce_col(cur);
        if (lead != -1) {
            u64 inv = M.inv(cur[0].second);
            for (auto& p : cur) p.second = M.mul(p.second, inv);
            has_piv[lead] = true;
            pivots[lead] = std::move(cur);
            rank_N++;
        }
    }

    // Test if A_src itself is already in the image of M_N!
    SparseCol colA;
    for (auto& [ma, ca] : poly_scale(poly_mul(A_src, D2, M), 2, M)) {
        colA.emplace_back(eq_idx.at(ma), ca);
    }
    norm_col(colA, M);
    bool hits_A = (reduce_col(colA) == -1);

    // Now insert all 40 odd-odd S columns!
    auto S_mons = gen_mons4(6, 1, 1);
    int rank_full = rank_N;
    for (const auto& ms : S_mons) {
        SparseCol cur;
        for (auto& [md, cd] : D2) {
            Mon4 mt{{ms.e[0]+md.e[0], ms.e[1]+md.e[1], ms.e[2]+md.e[2], ms.e[3]+md.e[3]}};
            cur.emplace_back(eq_idx.at(mt), M.mul(2, cd));
        }
        norm_col(cur, M);
        int lead = reduce_col(cur);
        if (lead != -1) {
            u64 inv = M.inv(cur[0].second);
            for (auto& p : cur) p.second = M.mul(p.second, inv);
            has_piv[lead] = true;
            pivots[lead] = std::move(cur);
            rank_full++;
        }
    }

    int dim_WD = 40 - rank_full + rank_N;
    return {degD, nRows, nN, rank_N, rank_full, dim_WD, dim_WD - base_E, hits_A};
}

static Poly4 make_H(u64 kappa, const ModCtx& M, bool coupled = true) {
    Poly4 H;
    H[Mon4{{6,0,0,0}}] = 1; H[Mon4{{0,6,0,0}}] = 1;
    H[Mon4{{0,0,6,0}}] = 1; H[Mon4{{0,0,0,6}}] = 1;
    H[Mon4{{5,1,0,0}}] = kappa % M.P;
    H[Mon4{{1,5,0,0}}] = M.mul(2, kappa);
    H[Mon4{{0,0,5,1}}] = kappa % M.P;
    H[Mon4{{0,0,1,5}}] = M.mul(3, kappa);
    if (coupled) {
        H[Mon4{{3,1,1,1}}] = kappa % M.P;
        H[Mon4{{1,1,3,1}}] = M.mul(2, kappa);
    }
    return H;
}

static int compute_poly_exact_dim(const Poly4& H, const ModCtx& M) {
    Poly4 dH[4];
    for (int i=0;i<4;++i) dH[i] = poly_deriv(H, i, M);
    auto S_mons = gen_mons4(6, 1, 1);
    std::map<Mon4, int> s_idx;
    for (int i=0;i<40;++i) s_idx[S_mons[i]] = i;

    std::vector<Poly4> gens;
    Poly4 v[4];
    for (int i=0;i<4;++i) { Mon4 m{{0,0,0,0}}; m.e[i]=1; v[i][m]=1; }
    for (int i : {0,1}) for (int j : {2,3}) gens.push_back(poly_mul(v[j], dH[i], M));
    for (int j : {2,3}) for (int i : {0,1}) gens.push_back(poly_mul(v[i], dH[j], M));

    std::vector<SparseCol> pivots(40);
    std::vector<bool> has_piv(40, false);
    int rk = 0;
    for (auto& g : gens) {
        SparseCol cur;
        for (auto& [m, c] : g) cur.emplace_back(s_idx.at(m), c);
        norm_col(cur, M);
        while (!cur.empty()) {
            int lead = cur[0].first;
            if (!has_piv[lead]) {
                u64 inv = M.inv(cur[0].second);
                for (auto& p : cur) p.second = M.mul(p.second, inv);
                has_piv[lead] = true; pivots[lead] = std::move(cur); rk++; break;
            }
            u64 fac = cur[0].second;
            SparseCol tmp;
            for (size_t a=1;a<cur.size();++a) tmp.push_back(cur[a]);
            for (size_t b=1;b<pivots[lead].size();++b)
                tmp.emplace_back(pivots[lead][b].first, M.sub(0, M.mul(fac, pivots[lead][b].second)));
            norm_col(tmp, M);
            cur.swap(tmp);
        }
    }
    return rk;
}

int main() {
    // 1. Rational covariants & baseline candidates at TWO independent large primes
    ModCtx M1{1000000007ULL}, M2{998244353ULL};
    u64 kap1 = M1.mul(5, M1.inv(7)), kap2 = M2.mul(5, M2.inv(7));
    Poly4 H1 = make_H(kap1, M1, true), H2 = make_H(kap2, M2, true);
    Poly4 A_src{{Mon4{{3,0,3,0}}, 1}};

    auto build_cands = [](const ModCtx& M, u64 kap) {
        Poly4 H = make_H(kap, M, true);
        Poly4 alpha, beta, gamma, delta;
        alpha[Mon4{{0,0,6,0}}] = 1; alpha[Mon4{{0,0,5,1}}] = kap;
        alpha[Mon4{{0,0,1,5}}] = M.mul(3, kap); alpha[Mon4{{0,0,0,6}}] = 1;
        beta[Mon4{{1,1,3,1}}] = M.mul(2, kap);
        gamma[Mon4{{3,1,1,1}}] = kap;
        delta[Mon4{{6,0,0,0}}] = 1; delta[Mon4{{5,1,0,0}}] = kap;
        delta[Mon4{{1,5,0,0}}] = M.mul(2, kap); delta[Mon4{{0,6,0,0}}] = 1;

        Poly4 h2 = poly_sub(poly_mul(beta, beta, M), poly_scale(poly_mul(alpha, gamma, M), 3, M), M);
        Poly4 h1 = poly_sub(poly_mul(beta, gamma, M), poly_scale(poly_mul(alpha, delta, M), 9, M), M);
        Poly4 h0 = poly_sub(poly_mul(gamma, gamma, M), poly_scale(poly_mul(beta, delta, M), 3, M), M);
        Poly4 H_cub = poly_add(h2, poly_add(h1, h0, M), M);

        Poly4 disc_cub = poly_mul(poly_mul(beta, beta, M), poly_mul(gamma, gamma, M), M);
        disc_cub = poly_sub(disc_cub, poly_scale(poly_mul(alpha, poly_pow(gamma, 3, M), M), 4, M), M);
        disc_cub = poly_sub(disc_cub, poly_scale(poly_mul(poly_pow(beta, 3, M), delta, M), 4, M), M);
        disc_cub = poly_sub(disc_cub, poly_scale(poly_mul(poly_mul(alpha, alpha, M), poly_mul(delta, delta, M), M), 27, M), M);
        disc_cub = poly_add(disc_cub, poly_scale(poly_mul(poly_mul(alpha, beta, M), poly_mul(gamma, delta, M), M), 18, M), M);

        Poly4 va{{Mon4{{1,0,0,0}}, 1}}, vb{{Mon4{{0,1,0,0}}, 1}}, vc{{Mon4{{0,0,1,0}}, 1}}, vd{{Mon4{{0,0,0,1}}, 1}};
        Poly4 a2_2c2 = poly_add(poly_pow(va, 2, M), poly_scale(poly_pow(vc, 2, M), 2, M), M);
        Poly4 abcd = poly_mul(poly_mul(va, vb, M), poly_mul(vc, vd, M), M);

        Poly4 Hmat[4][4];
        for (int i=0;i<4;++i) for (int j=0;j<4;++j) Hmat[i][j] = poly_deriv(poly_deriv(H, i, M), j, M);
        auto det3 = [&](int r0, int r1, int r2, int c0, int c1, int c2) {
            int rs[3] = {r0,r1,r2}, cs[3] = {c0,c1,c2};
            Poly4 t0 = poly_mul(Hmat[rs[0]][cs[0]], poly_sub(poly_mul(Hmat[rs[1]][cs[1]], Hmat[rs[2]][cs[2]], M), poly_mul(Hmat[rs[1]][cs[2]], Hmat[rs[2]][cs[1]], M), M), M);
            Poly4 t1 = poly_mul(Hmat[rs[0]][cs[1]], poly_sub(poly_mul(Hmat[rs[1]][cs[0]], Hmat[rs[2]][cs[2]], M), poly_mul(Hmat[rs[1]][cs[2]], Hmat[rs[2]][cs[0]], M), M), M);
            Poly4 t2 = poly_mul(Hmat[rs[0]][cs[2]], poly_sub(poly_mul(Hmat[rs[1]][cs[0]], Hmat[rs[2]][cs[1]], M), poly_mul(Hmat[rs[1]][cs[1]], Hmat[rs[2]][cs[0]], M), M), M);
            return poly_add(poly_sub(t0, t1, M), t2, M);
        };
        Poly4 Hess4 = poly_sub(
            poly_add(poly_mul(Hmat[0][0], det3(1,2,3, 1,2,3), M), poly_mul(Hmat[0][2], det3(1,2,3, 0,1,3), M), M),
            poly_add(poly_mul(Hmat[0][1], det3(1,2,3, 0,2,3), M), poly_mul(Hmat[0][3], det3(1,2,3, 0,1,2), M), M), M
        );

        Poly4 qsum = poly_add(poly_add(poly_pow(va,2,M), poly_pow(vb,2,M), M), poly_add(poly_pow(vc,2,M), poly_pow(vd,2,M), M), M);

        return std::vector<std::pair<std::string, Poly4>>{
            {"x0", va},
            {"x0*x1", poly_mul(va, vb, M)},
            {"x0^2+x1^2+y0^2+y1^2", qsum},
            {"a^2+2c^2", a2_2c2},
            {"abcd", abcd},
            {"abcd*(a^2+2c^2)", poly_mul(abcd, a2_2c2, M)},
            {"H", H},
            {"H^2", poly_pow(H, 2, M)},
            {"f_a", delta},
            {"f_c", alpha},
            {"f_a - f_c", poly_sub(delta, alpha, M)},
            {"f_a * f_c", poly_mul(delta, alpha, M)},
            {"H * f_a * f_c", poly_mul(H, poly_mul(delta, alpha, M), M)},
            {"h1 (deg 12)", h1},
            {"H_cub (deg 12)", H_cub},
            {"H * h1 (deg 18)", poly_mul(H, h1, M)},
            {"H * H_cub (deg 18)", poly_mul(H, H_cub, M)},
            {"H * a * c * H_cub (deg 20)", poly_mul(poly_mul(H, poly_mul(va, vc, M), M), H_cub, M)},
            {"Hess4(H) (deg 16)", Hess4},
            {"H * Hess4(H) (deg 22)", poly_mul(H, Hess4, M)},
            {"disc_cub (deg 24)", disc_cub},
            {"H_cub * f_a * f_c (deg 24)", poly_mul(H_cub, poly_mul(delta, alpha, M), M)},
            {"H * abcd * (a^2+2c^2) (deg 12)", poly_mul(H, poly_mul(abcd, a2_2c2, M), M)}
        };
    };

    auto c1 = build_cands(M1, kap1);
    auto c2 = build_cands(M2, kap2);
    std::cout << "=== PART 1: W_D / E AUDIT ACROSS RATIONAL COVARIANTS & STRUCTURED SUPPORTS ===\n";
    for (size_t i = 0; i < c1.size(); ++i) {
        auto r1 = audit_WD(c1[i].second, H1, A_src, M1, 8);
        auto r2 = audit_WD(c2[i].second, H2, A_src, M2, 8);
        std::cout << "  [" << c1[i].first << "] deg=" << r1.degD
                  << " | shape=" << r1.nRows << "x" << (40 + r1.nN)
                  << " | p1: dimW_D=" << r1.dim_WD << " (new=" << r1.new_classes << ")"
                  << " | p2: dimW_D=" << r2.dim_WD << " (new=" << r2.new_classes << ")\n";
    }

    // 2. POSITIVE CONTROL & NODAL-CONTACT TESTS AT p = 23 (where kappa = 5/7 = 4 mod 23!)
    ModCtx M23{23ULL};
    u64 kap23 = M23.mul(5, M23.inv(7)); // 4 mod 23
    Poly4 H23_uncoupled = make_H(kap23, M23, false);
    Poly4 H23_coupled   = make_H(kap23, M23, true);
    int E_unc = compute_poly_exact_dim(H23_uncoupled, M23);
    int E_cpl = compute_poly_exact_dim(H23_coupled, M23);
    std::cout << "\n=== PART 2: POSITIVE CONTROL & NODAL/ELLIPTIC CONTACT SUPPORTS AT p = 23 (kappa = " << kap23 << ") ===\n";
    std::cout << "  dim(E) mod 23: uncoupled H_0 -> " << E_unc << ", coupled H -> " << E_cpl << "\n";

    // Search for ALL binary contacts f_a(a,b) = C_a^3 - B_a^2 and f_c(c,d) = C_c^3 - B_c^2 mod 23!
    // Wait: f_c = C_c^3 - B_c^2 has known solution C_c = [15, 18, 16], B_c = [4, 6, 17, 22]!
    // Let's also check if f_a = C_a^3 - B_a^2 (or C_a^3 + R_a^2) has solutions mod 23 (or mod p where -1 is a square, e.g. p = 29 or 37 or 53)!
    Poly4 Cc, Bc;
    Cc[Mon4{{0,0,2,0}}] = 15; Cc[Mon4{{0,0,1,1}}] = 18; Cc[Mon4{{0,0,0,2}}] = 16;
    Bc[Mon4{{0,0,3,0}}] = 4;  Bc[Mon4{{0,0,2,1}}] = 6;  Bc[Mon4{{0,0,1,2}}] = 17; Bc[Mon4{{0,0,0,3}}] = 22;

    Poly4 va{{Mon4{{1,0,0,0}}, 1}}, vb{{Mon4{{0,1,0,0}}, 1}};
    Poly4 a3 = poly_pow(va, 3, M23), b3 = poly_pow(vb, 3, M23);
    Poly4 ab = poly_mul(va, vb, M23);

    std::vector<std::pair<std::string, Poly4>> nodal23 = {
        {"Cone a^3 + B_c(c,d) (deg 3)", poly_add(a3, Bc, M23)},
        {"Cone b^3 + B_c(c,d) (deg 3)", poly_add(b3, Bc, M23)},
        {"(a^3 + B_c)*(a^3 - B_c) = a^6 - B_c^2 (deg 6)", poly_sub(poly_pow(va,6,M23), poly_mul(Bc,Bc,M23), M23)},
        {"(a^3 + b^3)^2 - B_c^2 (deg 6)", poly_sub(poly_pow(poly_add(a3,b3,M23),2,M23), poly_mul(Bc,Bc,M23), M23)},
        {"Nodal contact H - C_c(c,d)^3 (deg 6)", poly_sub(H23_coupled, poly_pow(Cc, 3, M23), M23)},
        {"ab * (H - C_c^3) (deg 8)", poly_mul(ab, poly_sub(H23_coupled, poly_pow(Cc, 3, M23), M23), M23)}
    };

    std::cout << "  --- On UNCOUPLED H_0 = f_a + f_c (base E_0 = " << E_unc << "): ---\n";
    for (auto& [lbl, pol] : nodal23) {
        auto r = audit_WD(pol, H23_uncoupled, A_src, M23, E_unc);
        std::cout << "    [" << lbl << "] deg=" << r.degD << " -> dimW_D=" << r.dim_WD
                  << ", NEW_CLASSES=" << r.new_classes << ", hits_A=" << r.hits_A << "\n";
    }

    std::cout << "  --- On COUPLED H (base E = " << E_cpl << "): ---\n";
    for (auto& [lbl, pol] : nodal23) {
        auto r = audit_WD(pol, H23_coupled, A_src, M23, E_cpl);
        std::cout << "    [" << lbl << "] deg=" << r.degD << " -> dimW_D=" << r.dim_WD
                  << ", NEW_CLASSES=" << r.new_classes << ", hits_A=" << r.hits_A << "\n";
    }

    // Also find a prime p where BOTH f_a(a,b) = C_a^3 - B_a^2 AND f_c(c,d) = C_c^3 - B_c^2 have F_p solutions
    // so we can test benediamond's EXACT positive control D_bin = B_a(a,b)^2 - B_c(c,d)^2 on H_0 AND on H!
    for (u64 p : {29ULL, 31ULL, 37ULL, 41ULL, 43ULL, 47ULL, 53ULL, 59ULL, 61ULL, 67ULL, 71ULL, 73ULL, 79ULL, 83ULL, 89ULL, 97ULL, 101ULL, 103ULL, 107ULL, 109ULL, 113ULL}) {
        ModCtx Mp{p};
        u64 kp = Mp.mul(5, Mp.inv(7));
        // Brute-force C = c0 u^2 + c1 u v + c2 v^2 (3*p^2 -> p^3 = ~10^6? Wait: for p <= 53, p^3 <= 148877, instant!)
        if (p > 53) break;
        auto find_CB = [&](u64 k_lin5, u64 k_lin1) -> std::pair<std::vector<u64>, std::vector<u64>> {
            // Target binary sextic: u^6 + k_lin5 u^5 v + k_lin1 u v^5 + v^6 = C(u,v)^3 - B(u,v)^2
            for (u64 c0 = 0; c0 < p; ++c0)
            for (u64 c1 = 0; c1 < p; ++c1)
            for (u64 c2 = 0; c2 < p; ++c2) {
                // Compute G(t) = (c0 + c1 t + c2 t^2)^3 - (1 + k_lin5 t + k_lin1 t^5 + t^6) = B(t)^2
                u64 C[3] = {c0, c1, c2};
                u64 C2[5] = {0}, C3[7] = {0};
                for (int i=0;i<3;++i) for (int j=0;j<3;++j) C2[i+j] = Mp.add(C2[i+j], Mp.mul(C[i], C[j]));
                for (int i=0;i<5;++i) for (int j=0;j<3;++j) C3[i+j] = Mp.add(C3[i+j], Mp.mul(C2[i], C[j]));
                u64 G[7];
                G[0] = Mp.sub(C3[0], 1);
                G[1] = Mp.sub(C3[1], k_lin5);
                G[2] = C3[2];
                G[3] = C3[3];
                G[4] = C3[4];
                G[5] = Mp.sub(C3[5], k_lin1);
                G[6] = Mp.sub(C3[6], 1);
                // Check if G(t) = (b0 + b1 t + b2 t^2 + b3 t^3)^2 with G != 0 and gcd(C, B) = 1
                // Check b0^2 = G[0], b3^2 = G[6]
                for (u64 b0 = 0; b0 < p; ++b0) {
                    if (Mp.mul(b0, b0) != G[0]) continue;
                    if (b0 == 0 && G[1] != 0) continue;
                    if (b0 == 0) continue; // assume b0 != 0 for fast unique extraction
                    u64 inv2b0 = Mp.inv(Mp.mul(2, b0));
                    u64 b1 = Mp.mul(G[1], inv2b0);
                    u64 b2 = Mp.mul(Mp.sub(G[2], Mp.mul(b1, b1)), inv2b0);
                    u64 b3 = Mp.mul(Mp.sub(G[3], Mp.mul(2, Mp.mul(b1, b2))), inv2b0);
                    if (Mp.add(Mp.mul(b2, b2), Mp.mul(2, Mp.mul(b1, b3))) != G[4]) continue;
                    if (Mp.mul(2, Mp.mul(b2, b3)) != G[5]) continue;
                    if (Mp.mul(b3, b3) != G[6]) continue;
                    if (c0 == 0 && c1 == 0 && c2 == 0) continue;
                    return {{c0, c1, c2}, {b0, b1, b2, b3}};
                }
            }
            return {{}, {}};
        };
        auto [Ca, Ba] = find_CB(kp, Mp.mul(2, kp));
        auto [Cc_p, Bc_p] = find_CB(kp, Mp.mul(3, kp));
        if (!Ba.empty() && !Bc_p.empty()) {
            std::cout << "\n  *** FOUND SIMULTANEOUS BINARY CONTACTS AT p = " << p << " ***\n";
            Poly4 pBa, pBc, pCa, pCc;
            for (int j=0;j<4;++j) { pBa[Mon4{{3-j,j,0,0}}] = Ba[j]; pBc[Mon4{{0,0,3-j,j}}] = Bc_p[j]; }
            for (int j=0;j<3;++j) { pCa[Mon4{{2-j,j,0,0}}] = Ca[j]; pCc[Mon4{{0,0,2-j,j}}] = Cc_p[j]; }
            Poly4 D_bin = poly_sub(poly_mul(pBa, pBa, Mp), poly_mul(pBc, pBc, Mp), Mp);
            Poly4 H0_p = make_H(kp, Mp, false);
            Poly4 H_p  = make_H(kp, Mp, true);
            int e0 = compute_poly_exact_dim(H0_p, Mp);
            int e1 = compute_poly_exact_dim(H_p, Mp);
            auto r_unc = audit_WD(D_bin, H0_p, A_src, Mp, e0);
            auto r_cpl = audit_WD(D_bin, H_p,  A_src, Mp, e1);
            std::cout << "    D_bin = B_a^2 - B_c^2 (deg 6) on UNCOUPLED H_0: dimW_D=" << r_unc.dim_WD
                      << ", NEW_CLASSES=" << r_unc.new_classes << ", hits_A=" << r_unc.hits_A << "\n";
            std::cout << "    D_bin = B_a^2 - B_c^2 (deg 6) on COUPLED H    : dimW_D=" << r_cpl.dim_WD
                      << ", NEW_CLASSES=" << r_cpl.new_classes << ", hits_A=" << r_cpl.hits_A << "\n";
            Poly4 D_sumC = poly_add(pCa, pCc, Mp);
            auto r_sc0 = audit_WD(D_sumC, H0_p, A_src, Mp, e0);
            auto r_sc1 = audit_WD(D_sumC, H_p,  A_src, Mp, e1);
            std::cout << "    C_a + C_c (deg 2) on H_0: new=" << r_sc0.new_classes
                      << " | on H: new=" << r_sc1.new_classes << "\n";
        }
    }
    return 0;
}
