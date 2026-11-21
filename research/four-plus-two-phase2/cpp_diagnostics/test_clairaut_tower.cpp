// Exact Clairaut-Koszul Obstruction Tower for m = 0, deg(D) = 6!
// Constructs R^(2)(D) (quadratics), R^(3)(D) (cubics), R^(4)(D) (quartics), R^(5)(D) (quintics)
// in PURE D-VARIABLES ALONE (zero K-matrices!) and tests both:
//   (a) Complex Gauss-Newton on P^{nD-1} (PURE_(1,1), PURE_(0,0), and FULL_ALL_49!), and
//   (b) Exact F_p Macaulay rank in pure D^d monomials!

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <unordered_map>
#include <vector>

using u64 = uint64_t;
using u32 = uint32_t;
using cd = std::complex<double>;
static constexpr u64 MOD = 1000000007ULL;

static inline u32 mod_mul(u32 a, u32 b) { return (u32)((u64)a * b % MOD); }
static inline u32 mod_sub(u32 a, u32 b) { return a >= b ? a - b : (u32)(a + MOD - b); }
static inline u32 mod_add(u32 a, u32 b) { u32 s = a + b; return s >= MOD ? (u32)(s - MOD) : s; }
static u32 mod_pow(u32 a, u64 e) {
    u32 r = 1;
    while (e) { if (e & 1) r = mod_mul(r, a); a = mod_mul(a, a); e >>= 1; }
    return r;
}
static inline u32 mod_inv(u32 a) { return mod_pow(a, MOD - 2); }

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

using Poly4 = std::map<Mon4, u32>;
static void poly_add_term(Poly4& p, const Mon4& m, u32 c) {
    if (!c) return;
    auto it = p.find(m);
    if (it == p.end()) p[m] = c;
    else { u32 s = mod_add(it->second, c); if (!s) p.erase(it); else it->second = s; }
}
static Poly4 poly_mul(const Poly4& A, const Poly4& B) {
    Poly4 C;
    for (auto& [ma, ca] : A)
        for (auto& [mb, cb] : B)
            poly_add_term(C, Mon4{{ma.e[0]+mb.e[0], ma.e[1]+mb.e[1], ma.e[2]+mb.e[2], ma.e[3]+mb.e[3]}}, mod_mul(ca, cb));
    return C;
}
static Poly4 poly_deriv(const Poly4& A, int v) {
    Poly4 C;
    for (auto& [m, c] : A) if (m.e[v] > 0) {
        Mon4 md = m; u32 ev = (u32)md.e[v]; md.e[v]--;
        poly_add_term(C, md, mod_mul(c, ev));
    }
    return C;
}

// Build linear splitting of C[x]_12 = J(H)_12 (+) W_12 over F_MOD:
// Given any P in C[x]_12 (vector of length 455), returns:
//   - Pi(P) in W_12 (vector of length 35: non-pivot monomials of J(H)_12!)
//   - L_i(P) in C[x]_7 (4 polynomials of deg 7) such that P = sum_i L_i(P) d_i H + Pi(P)!
struct JH12Splitter {
    std::vector<Mon4> m12, m7;
    std::map<Mon4, int> idx12, idx7;
    std::vector<int> free_rows;       // 35 indices in 0..454 forming W_12
    // For each input monomial r in 0..454:
    //   if r is in free_rows: Pi(e_r) = e_r, L(e_r) = 0? NO! Wait:
    //   Let M (455 x 480) be multiplication (g_0..g_3) -> sum g_i d_i H!
    //   Column-reduce M^T (or row-reduce M with column pivots? Wait: we want M * L_mat = I - Pi_mat!):
    //   Do RREF on [M^T (480 x 455)]? Let's be super simple:
    //   Find 420 independent columns of M (out of 480) and 420 pivot rows of M!
    //   Let S_rows (size 420) be the pivot rows of M, and W_rows (size 35) be the OTHER 35 rows of M!
    //   Wait: if M_{S_rows, S_cols} (420 x 420) is invertible, does C[x]_12 = Im(M) (+) span(e_{W_rows})?!
    //   YES!! Because if sum_{w in W_rows} c_w e_w = M u, then looking at S_rows gives 0 = M_{S_rows, .} u;
    //   Wait: ONLY if rank(M) == 420! Since dim(R / J(H))_12 == 35 and 455 - 35 == 420, rank(M) IS 420,
    //   SO ANY 420 INDEPENDENT ROWS OF M CORRESPOND TO 35 COMPLEMENTARY BASIS MONOMIALS SPANNING R/J(H)_12!!
    //   Let's check: if u = M_{S_rows, S_cols}^{-1} P_{S_rows} (extended by 0 outside S_cols),
    //   then P - M u vanishes IDENTICALLY on all 420 S_rows, so P - M u is supported PURELY on W_rows (size 35)!!
    //   SO: L(P) := u, and Pi(P) := (P - M u)|_{W_rows}!!
    //   CHECK: M * L(P) + sum_{w in W_rows} Pi(P)_w e_w = M u + (P - M u) = P IDENTICALLY!!
    std::vector<int> piv_rows, w_rows, piv_cols;
    // L_op[r] = sparse list of (col_in_480, coeff) for r in 0..454
    // Pi_op[r] = list of (w_idx in 0..34, coeff) for r in 0..454
    std::vector<std::vector<std::pair<int, u32>>> L_op;
    std::vector<std::vector<std::pair<int, u32>>> Pi_op;

    void init(const Poly4 dH[4]) {
        m12 = gen_mons4(12);
        m7  = gen_mons4(7);
        for (int i=0;i<(int)m12.size();++i) idx12[m12[i]] = i;
        for (int i=0;i<(int)m7.size();++i)  idx7[m7[i]]   = i;
        int nr = (int)m12.size(); // 455
        int n7 = (int)m7.size();  // 120
        int nc = 4 * n7;          // 480

        std::vector<std::vector<u32>> M(nr, std::vector<u32>(nc, 0));
        for (int i=0;i<4;++i)
            for (int j=0;j<n7;++j)
                for (auto& [m, c] : poly_mul(Poly4{{m7[j], 1}}, dH[i]))
                    M[idx12.at(m)][i * n7 + j] = mod_add(M[idx12.at(m)][i * n7 + j], c);

        // First select 420 independent columns of M using RREF on M
        auto Mcopy = M;
        int r = 0;
        for (int c=0; c<nc && r<nr; ++c) {
            int sel=-1;
            for (int i=r;i<nr;++i) if (Mcopy[i][c]) { sel=i; break; }
            if (sel==-1) continue;
            std::swap(Mcopy[r], Mcopy[sel]);
            u32 inv = mod_inv(Mcopy[r][c]);
            for (int j=c;j<nc;++j) Mcopy[r][j] = mod_mul(Mcopy[r][j], inv);
            for (int i=r+1;i<nr;++i) if (Mcopy[i][c]) {
                u32 f = Mcopy[i][c];
                for (int j=c;j<nc;++j) Mcopy[i][j] = mod_sub(Mcopy[i][j], mod_mul(f, Mcopy[r][j]));
            }
            piv_cols.push_back(c);
            r++;
        }
        int rk = r; // should be 420!
        // Now select 420 independent ROWS of M[:, piv_cols] by row-reducing [M[:, piv_cols] | I_455]!
        std::vector<std::vector<u32>> Aug(nr, std::vector<u32>(rk + nr, 0));
        for (int i=0;i<nr;++i) {
            for (int j=0;j<rk;++j) Aug[i][j] = M[i][piv_cols[j]];
            Aug[i][rk + i] = 1;
        }
        std::vector<int> row_order(nr);
        for (int i=0;i<nr;++i) row_order[i] = i;
        int rr = 0;
        for (int c=0; c<rk && rr<nr; ++c) {
            int sel=-1;
            for (int i=rr;i<nr;++i) if (Aug[i][c]) { sel=i; break; }
            if (sel==-1) continue;
            std::swap(Aug[rr], Aug[sel]);
            std::swap(row_order[rr], row_order[sel]);
            u32 inv = mod_inv(Aug[rr][c]);
            for (int j=c;j<rk+nr;++j) Aug[rr][j] = mod_mul(Aug[rr][j], inv);
            for (int i=0;i<nr;++i) if (i!=rr && Aug[i][c]) {
                u32 f = Aug[i][c];
                for (int j=c;j<rk+nr;++j) Aug[i][j] = mod_sub(Aug[i][j], mod_mul(f, Aug[rr][j]));
            }
            rr++;
        }
        // Now Aug[0..rk-1, 0..rk-1] == I_rk, so Aug[0..rk-1, rk..rk+nr-1] gives u = M_{S,S}^{-1} P_S!
        // And Aug[rk..nr-1, 0..rk-1] == 0, while row_order[rk..nr-1] are the 35 w_rows!
        // Wait: even simpler! Look at Aug[rk + w, rk + r] for w = 0 .. nr-rk-1:
        // Since row operations zeroed columns 0..rk-1 of Aug[rk+w], Aug[rk+w, rk..rk+nr-1] * P IS LITERALLY (P - M u)_{w_rows[w]}!!
        for (int i=0;i<rk;++i) piv_rows.push_back(row_order[i]);
        for (int i=rk;i<nr;++i) w_rows.push_back(row_order[i]);

        L_op.resize(nr);
        Pi_op.resize(nr);
        for (int in_r = 0; in_r < nr; ++in_r) {
            for (int j = 0; j < rk; ++j) {
                u32 v = Aug[j][rk + in_r];
                if (v) L_op[in_r].emplace_back(piv_cols[j], v);
            }
            for (int w = 0; w < nr - rk; ++w) {
                u32 v = Aug[rk + w][rk + in_r];
                if (v) Pi_op[in_r].emplace_back(w, v);
            }
        }
        std::cout << "[JH12Splitter] nr=" << nr << ", rk(J(H)_12)=" << rk << ", dim(W_12)=" << w_rows.size() << "\n";
    }
};

static void gen_tuples_rec(int n, int d, int start, std::vector<int>& cur, std::vector<std::vector<int>>& out) {
    if ((int)cur.size() == d) { out.push_back(cur); return; }
    for (int i = start; i < n; ++i) {
        cur.push_back(i);
        gen_tuples_rec(n, d, i, cur, out);
        cur.pop_back();
    }
}
struct VectorHasher {
    size_t operator()(const std::vector<int>& v) const {
        size_t h = 1469598103934665603ULL;
        for (int x : v) { h ^= (size_t)(x + 1); h *= 1099511628211ULL; }
        return h;
    }
};

static int dense_reduce_inplace(std::vector<u64>& B, int nrows, int ncols) {
    int rd = 0;
    for (int c = 0; c < ncols && rd < nrows; ++c) {
        int piv = -1;
        for (int i = rd; i < nrows; ++i) {
            u64 v = B[(size_t)i * ncols + c] %= MOD;
            if (v != 0) { piv = i; break; }
        }
        if (piv == -1) continue;
        if (piv != rd) for (int j = c; j < ncols; ++j) std::swap(B[(size_t)rd * ncols + j], B[(size_t)piv * ncols + j]);
        u32 inv = mod_inv((u32)(B[(size_t)rd * ncols + c] % MOD));
        u64* pr = &B[(size_t)rd * ncols];
        for (int j = c; j < ncols; ++j) pr[j] = (pr[j] % MOD) * inv % MOD;
        for (int i = rd + 1; i < nrows; ++i) {
            u64 f = B[(size_t)i * ncols + c] % MOD;
            if (!f) { B[(size_t)i * ncols + c] = 0; continue; }
            u64 neg_f = MOD - f;
            u64* cr = &B[(size_t)i * ncols];
            cr[c] = 0;
            for (int j = c + 1; j < ncols; ++j) cr[j] += neg_f * pr[j];
        }
        rd++;
        if ((rd & 15) == 0) {
            for (int i = rd; i < nrows; ++i) {
                u64* cr = &B[(size_t)i * ncols];
                for (int j = c + 1; j < ncols; ++j) cr[j] %= MOD;
            }
        }
    }
    B.resize((size_t)rd * ncols);
    for (u64& z : B) z %= MOD;
    return rd;
}

// Pointwise evaluator over F_MOD for ANY D vector in E1: evaluates R^(2)(D), R^(3)(D), R^(4)(D), R^(5)(D), R^(6)(D)!
// WAIT A SECOND!! Look at that:
// Not only can we expand R^(m)(D) symbolically in D^m monomials, we can ALSO evaluate R^(2)(D)..R^(k)(D)
// symbolically up to d = 5 in PURE_(1,1) (nD = 21 -> tup[4] = 10,626, tup[5] = 53,130!)!!
static void run_clairaut_sector(const std::string& name, int par_ab, int par_cd, const Poly4& H7, const Poly4& A7, const JH12Splitter& split, int max_d) {
    auto t0 = std::chrono::steady_clock::now();
    std::cout << "\n======================================================\n";
    std::cout << "[Clairaut-Koszul Tower] Sector " << name << " (max_d=" << max_d << ")\n";

    Poly4 dH[4];
    for (int i=0;i<4;++i) dH[i] = poly_deriv(H7, i);
    auto D_mons = gen_mons4(6, par_ab, par_cd);
    std::vector<Mon4> N_mons[4];
    for (int i=0;i<4;++i) {
        if (par_ab == -1) N_mons[i] = gen_mons4(7, -1, -1);
        else if (i < 2)   N_mons[i] = gen_mons4(7, par_ab, par_cd ^ 1);
        else              N_mons[i] = gen_mons4(7, par_ab ^ 1, par_cd);
    }
    int nD_raw = (int)D_mons.size(), off[5] = {0};
    for (int i=0;i<4;++i) off[i+1] = off[i] + (int)N_mons[i].size();
    int nN_raw = off[4], tot1 = nN_raw + nD_raw;

    auto eq12 = (par_ab == -1) ? gen_mons4(12, -1, -1) : gen_mons4(12, par_ab ^ 1, par_cd ^ 1);
    std::map<Mon4, int> eq12_idx;
    for (int i=0;i<(int)eq12.size();++i) eq12_idx[eq12[i]] = i;

    std::vector<std::vector<u32>> M1(eq12.size(), std::vector<u32>(tot1, 0));
    for (int i=0;i<4;++i)
        for (int j=0;j<(int)N_mons[i].size();++j)
            for (auto& [m, c] : poly_mul(Poly4{{N_mons[i][j], 1}}, dH[i]))
                M1[eq12_idx.at(m)][off[i]+j] = mod_add(M1[eq12_idx.at(m)][off[i]+j], c);
    for (int j=0;j<nD_raw;++j)
        for (auto& [m, c] : poly_mul(Poly4{{D_mons[j], 1}}, A7))
            M1[eq12_idx.at(m)][nN_raw+j] = mod_sub(M1[eq12_idx.at(m)][nN_raw+j], c);

    int nrows = (int)M1.size(), r = 0;
    std::vector<int> piv_col, col_to_row(tot1, -1);
    for (int c=0;c<tot1 && r<nrows;++c) {
        int sel=-1;
        for (int i=r;i<nrows;++i) if (M1[i][c]) { sel=i; break; }
        if (sel==-1) continue;
        std::swap(M1[r], M1[sel]);
        u32 inv = mod_inv(M1[r][c]);
        for (int j=c;j<tot1;++j) M1[r][j] = mod_mul(M1[r][j], inv);
        for (int i=0;i<nrows;++i) if (i!=r && M1[i][c]) {
            u32 f = M1[i][c];
            for (int j=c;j<tot1;++j) M1[i][j] = mod_sub(M1[i][j], mod_mul(f, M1[r][j]));
        }
        piv_col.push_back(c); col_to_row[c] = r++;
    }
    std::vector<std::vector<u32>> E1_vecs;
    for (int c=0;c<tot1;++c) if (col_to_row[c]==-1) {
        std::vector<u32> v(tot1, 0); v[c] = 1;
        for (int pr=0;pr<r;++pr) if (M1[pr][c]) v[piv_col[pr]] = mod_sub(0, M1[pr][c]);
        bool dz = true;
        for (int j=0;j<nD_raw;++j) if (v[nN_raw+j]) { dz=false; break; }
        if (!dz) E1_vecs.push_back(std::move(v));
    }
    int nD = (int)E1_vecs.size();
    std::cout << "  nD = " << nD << "\n";

    std::vector<Poly4> D_b(nD), dD_b[4], divN0(nD), N0_b[4];
    for (int i=0;i<4;++i) { dD_b[i].resize(nD); N0_b[i].resize(nD); }
    for (int a=0;a<nD;++a) {
        for (int j=0;j<nD_raw;++j) if (E1_vecs[a][nN_raw+j]) D_b[a][D_mons[j]] = E1_vecs[a][nN_raw+j];
        for (int i=0;i<4;++i) {
            dD_b[i][a] = poly_deriv(D_b[a], i);
            for (int j=0;j<(int)N_mons[i].size();++j) if (E1_vecs[a][off[i]+j]) N0_b[i][a][N_mons[i][j]] = E1_vecs[a][off[i]+j];
            for (auto& [m, c] : poly_deriv(N0_b[i][a], i)) poly_add_term(divN0[a], m, c);
        }
    }

    // Precompute linear operator T_a : C[x]_7^4 (size 480) -> C[x]_12 (size 455) for each a = 0..nD-1:
    // T_a(g_0..g_3) := D_b[a] * sum_i d_i(g_i) - sum_i g_i * dD_b[i][a]!
    int n7 = (int)split.m7.size(); // 120
    int nc7 = 4 * n7;              // 480
    int nr12 = (int)split.m12.size(); // 455
    // Combine T_a with split.L_op and split.Pi_op!
    // Step_F[a][col7] -> sparse list of (out_col7, coeff)
    // Step_R[a][col7] -> sparse list of (w_idx in 0..34, coeff)
    std::vector<std::vector<std::vector<std::pair<int, u32>>>> Step_F(nD, std::vector<std::vector<std::pair<int, u32>>>(nc7));
    std::vector<std::vector<std::vector<std::pair<int, u32>>>> Step_R(nD, std::vector<std::vector<std::pair<int, u32>>>(nc7));

    for (int a = 0; a < nD; ++a) {
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < n7; ++j) {
                int col7 = i * n7 + j;
                Poly4 g_ij{{split.m7[j], 1}};
                Poly4 dg = poly_deriv(g_ij, i);
                Poly4 P = poly_mul(D_b[a], dg);
                for (auto& [m, c] : poly_mul(g_ij, dD_b[i][a])) poly_add_term(P, m, mod_sub(0, c));
                // Now apply split.L_op and split.Pi_op to P!
                std::map<int, u32> f_acc, r_acc;
                for (auto& [m, c] : P) {
                    int r12 = split.idx12.at(m);
                    for (auto& [c_out, lval] : split.L_op[r12])
                        f_acc[c_out] = mod_add(f_acc[c_out], mod_mul(c, lval));
                    for (auto& [w_out, pval] : split.Pi_op[r12])
                        r_acc[w_out] = mod_add(r_acc[w_out], mod_mul(c, pval));
                }
                for (auto& [k, v] : f_acc) if (v) Step_F[a][col7].emplace_back(k, v);
                for (auto& [k, v] : r_acc) if (v) Step_R[a][col7].emplace_back(k, v);
            }
        }
    }

    // Tuple tables for degrees 2 .. max_d
    std::vector<std::vector<std::vector<int>>> tup(max_d + 1);
    std::vector<std::unordered_map<std::vector<int>, int, VectorHasher>> idx(max_d + 1);
    for (int d = 2; d <= max_d; ++d) {
        std::vector<int> cv;
        gen_tuples_rec(nD, d, 0, cv, tup[d]);
        idx[d].reserve(tup[d].size() * 2 + 1);
        for (int i=0;i<(int)tup[d].size();++i) idx[d][tup[d][i]] = i;
    }

    // Base d = 2: compute F^(2)[mon2][col7] and R^(2)[w][mon2]!
    int nM2 = (int)tup[2].size();
    int nW = (int)split.w_rows.size(); // 35
    // Store F_cur as vector of size nM_d, each entry is sparse/dense vector over col7 (0..479)
    std::vector<std::vector<u32>> F_cur(nM2, std::vector<u32>(nc7, 0));
    std::vector<u64> ideal_basis((size_t)nW * nM2, 0);

    for (int m2 = 0; m2 < nM2; ++m2) {
        int a1 = tup[2][m2][0], a2 = tup[2][m2][1];
        Poly4 P = poly_mul(D_b[a1], divN0[a2]);
        for (int i=0;i<4;++i) for (auto& [m,c] : poly_mul(N0_b[i][a2], dD_b[i][a1])) poly_add_term(P, m, mod_sub(0, c));
        if (a1 < a2) {
            for (auto& [m,c] : poly_mul(D_b[a2], divN0[a1])) poly_add_term(P, m, c);
            for (int i=0;i<4;++i) for (auto& [m,c] : poly_mul(N0_b[i][a1], dD_b[i][a2])) poly_add_term(P, m, mod_sub(0, c));
        }
        for (auto& [m, c] : P) {
            int r12 = split.idx12.at(m);
            for (auto& [c_out, lval] : split.L_op[r12])
                F_cur[m2][c_out] = mod_add(F_cur[m2][c_out], mod_mul(c, lval));
            for (auto& [w_out, pval] : split.Pi_op[r12])
                ideal_basis[(size_t)w_out * nM2 + m2] = (ideal_basis[(size_t)w_out * nM2 + m2] + (u64)c * pval) % MOD;
        }
    }
    int rk_cur = dense_reduce_inplace(ideal_basis, nW, nM2);
    std::cout << "  [d=2] R^(2) rank = " << rk_cur << " / " << nM2 << "\n";

    // Iterate d = 3 .. max_d!
    for (int d = 3; d <= max_d; ++d) {
        auto td = std::chrono::steady_clock::now();
        int nPrev = (int)tup[d-1].size();
        int nCur  = (int)tup[d].size();
        std::vector<int> trans((size_t)nPrev * nD);
        for (int mp = 0; mp < nPrev; ++mp) {
            for (int a = 0; a < nD; ++a) {
                auto k = tup[d-1][mp]; k.push_back(a); std::sort(k.begin(), k.end());
                trans[(size_t)mp * nD + a] = idx[d].at(k);
            }
        }

        // 1. Compute F_next (size nCur x nc7) and R_new (size nW x nCur) from F_cur!
        std::vector<std::vector<u64>> F_next64(nCur, std::vector<u64>(nc7, 0));
        std::vector<u64> R_new((size_t)nW * nCur, 0);

        for (int mp = 0; mp < nPrev; ++mp) {
            for (int col7 = 0; col7 < nc7; ++col7) {
                u32 coeff = F_cur[mp][col7];
                if (!coeff) continue;
                for (int a = 0; a < nD; ++a) {
                    int mc = trans[(size_t)mp * nD + a];
                    for (auto& [c_out, fv] : Step_F[a][col7])
                        F_next64[mc][c_out] += (u64)coeff * fv;
                    for (auto& [w_out, rv] : Step_R[a][col7])
                        R_new[(size_t)w_out * nCur + mc] += (u64)coeff * rv;
                }
            }
        }
        F_cur.assign(nCur, std::vector<u32>(nc7, 0));
        for (int mc = 0; mc < nCur; ++mc)
            for (int c7 = 0; c7 < nc7; ++c7)
                F_cur[mc][c7] = (u32)(F_next64[mc][c7] % MOD);
        std::vector<std::vector<u64>>().swap(F_next64);

        // 2. Lift previous ideal_basis (rk_cur rows of length nPrev) by multiplying by all nD variables + append R_new (nW rows)!
        int tot_rows = rk_cur * nD + nW;
        std::vector<u64> B_next((size_t)tot_rows * nCur, 0);
        for (int r_i = 0; r_i < rk_cur; ++r_i) {
            for (int mp = 0; mp < nPrev; ++mp) {
                u64 val = ideal_basis[(size_t)r_i * nPrev + mp];
                if (!val) continue;
                for (int a = 0; a < nD; ++a) {
                    B_next[(size_t)(r_i * nD + a) * nCur + trans[(size_t)mp * nD + a]] += val;
                }
            }
        }
        for (int w = 0; w < nW; ++w) {
            for (int mc = 0; mc < nCur; ++mc) {
                B_next[(size_t)(rk_cur * nD + w) * nCur + mc] = R_new[(size_t)w * nCur + mc] % MOD;
            }
        }
        ideal_basis.swap(B_next);
        rk_cur = dense_reduce_inplace(ideal_basis, tot_rows, nCur);
        double sd = std::chrono::duration<double>(std::chrono::steady_clock::now() - td).count();
        std::cout << "  [d=" << d << " (" << sd << "s)] rows=" << tot_rows
                  << " -> ideal rank = " << rk_cur << " / " << nCur
                  << " (nullity = " << (nCur - rk_cur) << ")\n";
        if (rk_cur == nCur) {
            std::cout << "  *** UNCONDITIONAL EXCLUSION PROVED AT DEGREE d = " << d << " FOR SECTOR " << name << "!! ***\n";
            break;
        }
    }
    std::cout << "Finished " << name << " in " << std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count() << "s.\n";
}

int main() {
    Poly4 H7, A7, dH[4];
    H7[Mon4{{6,0,0,0}}] = 7; H7[Mon4{{0,6,0,0}}] = 7;
    H7[Mon4{{0,0,6,0}}] = 7; H7[Mon4{{0,0,0,6}}] = 7;
    H7[Mon4{{5,1,0,0}}] = 5; H7[Mon4{{1,5,0,0}}] = 10;
    H7[Mon4{{0,0,5,1}}] = 5; H7[Mon4{{0,0,1,5}}] = 15;
    H7[Mon4{{3,1,1,1}}] = 5; H7[Mon4{{1,1,3,1}}] = 10;
    A7[Mon4{{3,0,3,0}}] = 7;
    for (int i=0;i<4;++i) dH[i] = poly_deriv(H7, i);

    JH12Splitter split;
    split.init(dH);

    run_clairaut_sector("PURE_(1,1)", 1, 1, H7, A7, split, 4);
    run_clairaut_sector("PURE_(0,0)", 0, 0, H7, A7, split, 4);
    return 0;
}
