// High-speed C++ solver for:
// 1. m = 1 (D = H * D0) exact Bidegree-(k,1) Macaulay test for general D0 of degrees d0 = 1, 2, 3, 4, 5!
//    Equation (after Euler gauge M . grad(H) = 0):
//      3 * (D0 * div(M) - M . grad(D0)) + 2 * A * D0^2 = 0,
//    where M_k = sum_j Omega^{jk} d_j H (deg(Omega) = d0 + 2) and deg(eq) = 2*d0 + 6.
//    Wait: ALSO checks general CoreCertificate(H * D0, N) WITHOUT assuming H | (N.grad(H))?
//    Wait: earlier we proved ANY CoreCertificate(H * D0, N) satisfies H | (N . grad(H))! Let's verify BOTH!
// 2. Fixed structured & covariant denominators D (up to degree 36+):
//    Tests linear solvability of CoreCertificate(D, N) over F_p (p = 1000000007).

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

using u64 = uint64_t;
using i64 = int64_t;
static constexpr u64 MOD = 1000000007ULL;

static inline u64 mod_add(u64 a, u64 b) { u64 s = a + b; return s >= MOD ? s - MOD : s; }
static inline u64 mod_sub(u64 a, u64 b) { return a >= b ? a - b : a + MOD - b; }
static inline u64 mod_mul(u64 a, u64 b) { return (unsigned __int128)a * b % MOD; }
static u64 mod_pow(u64 a, u64 e) {
    u64 r = 1;
    while (e) { if (e & 1) r = mod_mul(r, a); a = mod_mul(a, a); e >>= 1; }
    return r;
}
static inline u64 mod_inv(u64 a) { return mod_pow(a, MOD - 2); }

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

static void poly_add_term(Poly4& p, const Mon4& m, u64 c) {
    if (!c) return;
    auto it = p.find(m);
    if (it == p.end()) p[m] = c;
    else {
        u64 s = mod_add(it->second, c);
        if (s == 0) p.erase(it);
        else it->second = s;
    }
}

static Poly4 poly_add(const Poly4& A, const Poly4& B) {
    Poly4 C = A;
    for (auto& [m, c] : B) poly_add_term(C, m, c);
    return C;
}

static Poly4 poly_sub(const Poly4& A, const Poly4& B) {
    Poly4 C = A;
    for (auto& [m, c] : B) poly_add_term(C, m, mod_sub(0, c));
    return C;
}

static Poly4 poly_scale(const Poly4& A, u64 s) {
    Poly4 C;
    if (!s) return C;
    for (auto& [m, c] : A) poly_add_term(C, m, mod_mul(c, s));
    return C;
}

static Poly4 poly_mul(const Poly4& A, const Poly4& B) {
    Poly4 C;
    for (auto& [ma, ca] : A)
        for (auto& [mb, cb] : B) {
            Mon4 mc{{ma.e[0]+mb.e[0], ma.e[1]+mb.e[1], ma.e[2]+mb.e[2], ma.e[3]+mb.e[3]}};
            poly_add_term(C, mc, mod_mul(ca, cb));
        }
    return C;
}

static Poly4 poly_pow(Poly4 A, int e) {
    Poly4 R{{Mon4{{0,0,0,0}}, 1}};
    while (e > 0) {
        if (e & 1) R = poly_mul(R, A);
        e >>= 1;
        if (e) A = poly_mul(A, A);
    }
    return R;
}

static Poly4 poly_deriv(const Poly4& A, int var) {
    Poly4 C;
    for (auto& [m, c] : A) {
        if (m.e[var] > 0) {
            Mon4 md = m;
            u64 coeff = mod_mul(c, (u64)md.e[var]);
            md.e[var]--;
            poly_add_term(C, md, coeff);
        }
    }
    return C;
}

using SparseRow = std::vector<std::pair<int, u64>>;

static void sparse_normalize(SparseRow& row) {
    if (row.empty()) return;
    std::sort(row.begin(), row.end(), [](const auto& a, const auto& b){ return a.first < b.first; });
    size_t w = 0;
    for (size_t i = 0; i < row.size(); ++i) {
        if (row[i].second == 0) continue;
        if (w > 0 && row[w-1].first == row[i].first) {
            row[w-1].second = mod_add(row[w-1].second, row[i].second);
            if (row[w-1].second == 0) --w;
        } else {
            row[w++] = row[i];
        }
    }
    row.resize(w);
}

struct SparseElimStats {
    int rank_K = 0;
    int rank_D = 0;
};

static SparseElimStats sparse_elim_two_block(std::vector<SparseRow>& rows, int nK_cols, int nD_cols) {
    int total_cols = nK_cols + nD_cols;
    std::vector<SparseRow> pivots(total_cols);
    std::vector<bool> has_pivot(total_cols, false);
    SparseElimStats st{0, 0};

    std::vector<int> order(rows.size());
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](int a, int b){ return rows[a].size() < rows[b].size(); });

    std::vector<std::pair<int, u64>> tmp;
    for (int idx : order) {
        SparseRow cur = std::move(rows[idx]);
        while (!cur.empty()) {
            int lead = cur[0].first;
            if (!has_pivot[lead]) {
                u64 inv = mod_inv(cur[0].second);
                for (auto& p : cur) p.second = mod_mul(p.second, inv);
                has_pivot[lead] = true;
                pivots[lead] = std::move(cur);
                if (lead < nK_cols) st.rank_K++;
                else st.rank_D++;
                break;
            }
            u64 factor = cur[0].second;
            const SparseRow& piv = pivots[lead];
            tmp.clear();
            tmp.reserve(cur.size() + piv.size());
            size_t i = 1, j = 1;
            while (i < cur.size() && j < piv.size()) {
                if (cur[i].first < piv[j].first) tmp.push_back(cur[i++]);
                else if (piv[j].first < cur[i].first) {
                    u64 val = mod_sub(0, mod_mul(factor, piv[j].second));
                    if (val) tmp.emplace_back(piv[j].first, val);
                    ++j;
                } else {
                    u64 val = mod_sub(cur[i].second, mod_mul(factor, piv[j].second));
                    if (val) tmp.emplace_back(cur[i].first, val);
                    ++i; ++j;
                }
            }
            while (i < cur.size()) tmp.push_back(cur[i++]);
            while (j < piv.size()) {
                u64 val = mod_sub(0, mod_mul(factor, piv[j].second));
                if (val) tmp.emplace_back(piv[j].first, val);
                ++j;
            }
            cur.swap(tmp);
        }
        if (st.rank_D == nD_cols) break;
    }
    return st;
}

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

// Test m = 1 (D = H * D0) for general D0 of degree d0 in pure parity (par_ab, par_cd) or full (-1, -1)!
// Equation: 3 * (D0 * div(M) - M . grad(D0)) + 2 * A7 * D0^2 = 0, where M . grad(H7) = 0, deg(M) = d0 + 7.
// Note: M_k = sum_{j} Omega^{jk} d_j H7 with skew-symmetric Omega^{jk} of degree d0 + 2!
static void test_m1_macaulay(int d0, int par_ab, int par_cd, const Poly4& H7, const Poly4& A7, const std::vector<int>& ks) {
    std::cout << "\n======================================================\n";
    std::cout << "[m=1] Testing deg(D0) = " << d0 << " (deg(D) = " << (d0+6) << "), parity = (" << par_ab << "," << par_cd << ")\n";
    Poly4 dH[4];
    for (int i = 0; i < 4; ++i) dH[i] = poly_deriv(H7, i);

    std::vector<Mon4> D0_mons = gen_mons4(d0, par_ab, par_cd);
    int nD = (int)D0_mons.size();
    std::vector<Poly4> D0_basis(nD), dD0_basis[4];
    for (int i = 0; i < 4; ++i) dD0_basis[i].resize(nD);
    for (int a = 0; a < nD; ++a) {
        D0_basis[a][D0_mons[a]] = 1;
        for (int i = 0; i < 4; ++i) dD0_basis[i][a] = poly_deriv(D0_basis[a], i);
    }

    // Build Koszul syzygy basis M^{(b)} of degree d0 + 7 using skew pairs (j < k) times monomials of deg d0 + 2!
    // Check parity: we need D0 * div(M) - M . grad(D0) to match parity of A7 * D0^2, which is ALWAYS (1, 1) when D0 has pure parity!
    // Wait: if D0 has parity (pa, pc), then D0^2 has parity (0,0), so A7*D0^2 has parity (1,1).
    // Thus div(M) needs parity (pa ^ 1, pc ^ 1).
    // Since div(M) has parity of d_k M_k, M_{0,1} needs parity (pa, pc ^ 1), and M_{2,3} needs parity (pa ^ 1, pc)!
    // Let's check each skew pair (j < k) for M = m * (d_j H e_k - d_k H e_j):
    // - (0, 1): e_0, e_1 get d_1 H, -d_0 H (parity (1,0)). For M_{0,1} to have (pa, pc^1), m needs (pa ^ 1, pc ^ 1).
    // - (2, 3): e_2, e_3 get d_3 H, -d_2 H (parity (0,1)). For M_{2,3} to have (pa^1, pc), m needs (pa ^ 1, pc ^ 1).
    // - (0, 2), (0, 3), (1, 2), (1, 3): e_k (in 2,3) gets d_j H (parity (1,0)) -> needs (pa^1, pc), so m needs (pa, pc)!
    //   Check e_j (in 0,1): gets -d_k H (parity (0,1)) * m (pa, pc) = (pa, pc^1)! MATCHES!
    std::vector<Poly4> K_basis[4], divK;
    int pairs[6][2] = {{0,1},{2,3},{0,2},{0,3},{1,2},{1,3}};
    for (int p = 0; p < 6; ++p) {
        int j = pairs[p][0], k = pairs[p][1];
        int req_pa = -1, req_pc = -1;
        if (par_ab != -1) {
            if (p < 2) { req_pa = par_ab ^ 1; req_pc = par_cd ^ 1; }
            else       { req_pa = par_ab;     req_pc = par_cd; }
        }
        auto om_mons = gen_mons4(d0 + 2, req_pa, req_pc);
        for (const auto& m : om_mons) {
            Poly4 pm{{m, 1}};
            Poly4 mj = poly_sub(Poly4{}, poly_mul(pm, dH[k])); // -m * d_k H in slot j
            Poly4 mk = poly_mul(pm, dH[j]);                    // +m * d_j H in slot k
            for (int i = 0; i < 4; ++i) {
                if (i == j) K_basis[i].push_back(mj);
                else if (i == k) K_basis[i].push_back(mk);
                else K_basis[i].push_back(Poly4{});
            }
            Poly4 dv = poly_add(poly_deriv(mj, j), poly_deriv(mk, k));
            divK.push_back(std::move(dv));
        }
    }
    int nK = (int)divK.size();

    // Target equations of degree 2*d0 + 6
    auto eq_mons = (par_ab == -1) ? gen_mons4(2 * d0 + 6, -1, -1) : gen_mons4(2 * d0 + 6, 1, 1);
    std::map<Mon4, int> eq_idx;
    for (int i = 0; i < (int)eq_mons.size(); ++i) eq_idx[eq_mons[i]] = i;
    int m_eq = (int)eq_mons.size();
    std::cout << "  nD0 = " << nD << ", raw Kosz generators nK = " << nK << ", m_eq = " << m_eq << "\n";

    struct QTerm { int a1, a2; u64 val; };
    struct BTerm { int a, b; u64 val; };
    std::vector<std::vector<QTerm>> Q_eq(m_eq);
    std::vector<std::vector<BTerm>> B_eq(m_eq);

    // Quadratic part: 2 * A7 * D0_{a1} * D0_{a2} (with factor 2 more if a1 < a2)
    for (int a1 = 0; a1 < nD; ++a1) {
        for (int a2 = a1; a2 < nD; ++a2) {
            u64 mult = (a1 == a2) ? 2 : 4;
            Poly4 P = poly_scale(poly_mul(A7, poly_mul(D0_basis[a1], D0_basis[a2])), mult);
            for (auto& [m, c] : P) Q_eq[eq_idx.at(m)].push_back({a1, a2, c});
        }
    }
    // Bilinear part: 3 * (D0_a * divK_b - K_b . grad(D0_a))
    for (int a = 0; a < nD; ++a) {
        for (int b = 0; b < nK; ++b) {
            Poly4 P = poly_mul(D0_basis[a], divK[b]);
            for (int i = 0; i < 4; ++i)
                P = poly_sub(P, poly_mul(K_basis[i][b], dD0_basis[i][a]));
            P = poly_scale(P, 3);
            for (auto& [m, c] : P) B_eq[eq_idx.at(m)].push_back({a, b, c});
        }
    }

    int active_eq = 0;
    for (int r = 0; r < m_eq; ++r) if (!Q_eq[r].empty() || !B_eq[r].empty()) active_eq++;

    for (int k : ks) {
        auto tk0 = std::chrono::steady_clock::now();
        std::vector<std::vector<int>> tup_km2, tup_km1, tup_k, cur;
        std::vector<int> cvec;
        gen_tuples_rec(nD, k - 2, 0, cvec, tup_km2);
        gen_tuples_rec(nD, k - 1, 0, cvec, tup_km1);
        gen_tuples_rec(nD, k,     0, cvec, tup_k);

        std::unordered_map<std::vector<int>, int, VectorHasher> idx_km1, idx_k;
        for (int i = 0; i < (int)tup_km1.size(); ++i) idx_km1[tup_km1[i]] = i;
        for (int i = 0; i < (int)tup_k.size(); ++i)   idx_k[tup_k[i]] = i;

        int nK_cols = (int)tup_km1.size() * nK;
        int nD_cols = (int)tup_k.size();
        long long n_rows = (long long)active_eq * (long long)tup_km2.size();
        if (n_rows > 15000 || (nK_cols + nD_cols) > 15000) {
            std::cout << "  [Skipping k=" << k << ": rows=" << n_rows << ", cols=" << (nK_cols+nD_cols) << " too large]\n";
            break;
        }
        std::vector<SparseRow> mac_rows;
        mac_rows.reserve((size_t)n_rows);
        std::vector<int> key1, key2;
        for (const auto& mult : tup_km2) {
            std::vector<int> ins1(nD);
            for (int a = 0; a < nD; ++a) {
                key1 = mult; key1.push_back(a); std::sort(key1.begin(), key1.end());
                ins1[a] = idx_km1.at(key1);
            }
            for (int r = 0; r < m_eq; ++r) {
                if (Q_eq[r].empty() && B_eq[r].empty()) continue;
                SparseRow row;
                row.reserve(Q_eq[r].size() + B_eq[r].size());
                for (const auto& bt : B_eq[r]) row.emplace_back(ins1[bt.a] * nK + bt.b, bt.val);
                for (const auto& qt : Q_eq[r]) {
                    key2 = mult; key2.push_back(qt.a1); key2.push_back(qt.a2); std::sort(key2.begin(), key2.end());
                    row.emplace_back(nK_cols + idx_k.at(key2), qt.val);
                }
                sparse_normalize(row);
                if (!row.empty()) mac_rows.push_back(std::move(row));
            }
        }
        auto st = sparse_elim_two_block(mac_rows, nK_cols, nD_cols);
        double secs = std::chrono::duration<double>(std::chrono::steady_clock::now() - tk0).count();
        std::cout << "  k=" << k << " (" << secs << "s): rows=" << n_rows
                  << ", rank_K=" << st.rank_K << "/" << nK_cols
                  << ", rank_D=" << st.rank_D << "/" << nD_cols
                  << " (D^" << k << " nullity=" << (nD_cols - st.rank_D) << ")\n";
        if (st.rank_D == nD_cols) {
            std::cout << "  >>> EXCLUDED: NO NONZERO SOLUTION FOR m=1, d0=" << d0 << " IN THIS SECTOR! <<<\n";
            break;
        }
    }
}

// Test whether ANY N solves CoreCertificate(D, N) for a FIXED homogeneous polynomial D!
// Uses pure parity if D has pure parity, otherwise full space!
// Solves [N_cols | rhs_col] via sparse elimination: if rhs_col gets a pivot, INCONSISTENT (no solution);
// if rhs_col reduces to 0, SOLUTION FOUND!!
static bool test_fixed_D(const std::string& label, const Poly4& D, const Poly4& H7, const Poly4& A7) {
    if (D.empty()) return false;
    int degD = D.begin()->first.deg();
    int pa = D.begin()->first.pa();
    int pc = D.begin()->first.pc();
    bool pure = true;
    for (auto& [m, c] : D) {
        if (m.pa() != pa || m.pc() != pc) { pure = false; break; }
    }

    Poly4 dH[4], dD[4];
    for (int i = 0; i < 4; ++i) {
        dH[i] = poly_deriv(H7, i);
        dD[i] = poly_deriv(D, i);
    }

    std::vector<Mon4> N_mons[4];
    for (int i = 0; i < 4; ++i) {
        if (!pure) N_mons[i] = gen_mons4(degD + 1, -1, -1);
        else if (i < 2) N_mons[i] = gen_mons4(degD + 1, pa, pc ^ 1);
        else            N_mons[i] = gen_mons4(degD + 1, pa ^ 1, pc);
    }
    int off[5] = {0};
    for (int i = 0; i < 4; ++i) off[i+1] = off[i] + (int)N_mons[i].size();
    int nN = off[4];

    auto eq_mons = (!pure) ? gen_mons4(2 * degD + 6, -1, -1) : gen_mons4(2 * degD + 6, 1, 1);
    std::map<Mon4, int> eq_idx;
    for (int i = 0; i < (int)eq_mons.size(); ++i) eq_idx[eq_mons[i]] = i;
    int m_eq = (int)eq_mons.size();

    std::vector<SparseRow> rows(m_eq);

    // Precompute -2*D*dH[i] - 3*H7*dD[i] (multiplier of N_i) and 3*H7*D (multiplier of d_i N_i)
    Poly4 HD3 = poly_scale(poly_mul(H7, D), 3);
    Poly4 coeff_Ni[4];
    for (int i = 0; i < 4; ++i) {
        Poly4 t1 = poly_scale(poly_mul(D, dH[i]), 2);
        Poly4 t2 = poly_scale(poly_mul(H7, dD[i]), 3);
        coeff_Ni[i] = poly_sub(Poly4{}, poly_add(t1, t2));
    }

    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < (int)N_mons[i].size(); ++j) {
            int col = off[i] + j;
            Mon4 mn = N_mons[i][j];
            // Term 1: coeff_Ni[i] * mn
            for (auto& [mc, cc] : coeff_Ni[i]) {
                Mon4 mt{{mc.e[0]+mn.e[0], mc.e[1]+mn.e[1], mc.e[2]+mn.e[2], mc.e[3]+mn.e[3]}};
                rows[eq_idx.at(mt)].emplace_back(col, cc);
            }
            // Term 2: HD3 * d_i(mn)
            if (mn.e[i] > 0) {
                Mon4 dm = mn;
                u64 e_val = (u64)dm.e[i];
                dm.e[i]--;
                for (auto& [mc, cc] : HD3) {
                    Mon4 mt{{mc.e[0]+dm.e[0], mc.e[1]+dm.e[1], mc.e[2]+dm.e[2], mc.e[3]+dm.e[3]}};
                    rows[eq_idx.at(mt)].emplace_back(col, mod_mul(cc, e_val));
                }
            }
        }
    }

    // RHS term: + 2 * A7 * D^2 in column nN!
    Poly4 rhs = poly_scale(poly_mul(A7, poly_mul(D, D)), 2);
    for (auto& [m, c] : rhs) {
        rows[eq_idx.at(m)].emplace_back(nN, c);
    }

    std::vector<SparseRow> clean_rows;
    clean_rows.reserve(m_eq);
    for (auto& r : rows) {
        sparse_normalize(r);
        if (!r.empty()) clean_rows.push_back(std::move(r));
    }

    auto st = sparse_elim_two_block(clean_rows, nN, 1);
    bool ok = (st.rank_D == 0);
    std::cout << "  Fixed D [" << label << "] (deg=" << degD << ", terms=" << D.size()
              << ", nN=" << nN << ", eqs=" << m_eq << "): rank_N=" << st.rank_K
              << ", inconsistent=" << st.rank_D << (ok ? "  *** SOLUTION FOUND!! ***" : "") << "\n";
    return ok;
}

int main() {
    Poly4 H7, A7;
    H7[Mon4{{6,0,0,0}}] = 7; H7[Mon4{{0,6,0,0}}] = 7;
    H7[Mon4{{0,0,6,0}}] = 7; H7[Mon4{{0,0,0,6}}] = 7;
    H7[Mon4{{5,1,0,0}}] = 5; H7[Mon4{{1,5,0,0}}] = 10;
    H7[Mon4{{0,0,5,1}}] = 5; H7[Mon4{{0,0,1,5}}] = 15;
    H7[Mon4{{3,1,1,1}}] = 5; H7[Mon4{{1,1,3,1}}] = 10;
    A7[Mon4{{3,0,3,0}}] = 7;

    // Part 1: Test m = 1 Macaulay for d0 = 1, 2, 3, 4!
    test_m1_macaulay(1, -1, -1, H7, A7, {2, 3, 4});
    test_m1_macaulay(2, 0, 0, H7, A7, {2, 3, 4});
    test_m1_macaulay(2, 1, 1, H7, A7, {2, 3, 4});
    test_m1_macaulay(2, -1, -1, H7, A7, {2, 3, 4});
    test_m1_macaulay(3, 1, 0, H7, A7, {2, 3, 4});
    test_m1_macaulay(3, 0, 1, H7, A7, {2, 3, 4});
    test_m1_macaulay(4, 0, 0, H7, A7, {2, 3, 4});
    test_m1_macaulay(4, 1, 1, H7, A7, {2, 3, 4});

    // Part 2: Build cubic covariants of H7 viewed as cubic in lambda^2:
    // alpha = f_c(c,d), beta = 10 a b c^3 d, gamma = 5 a^3 b c d, delta = f_a(a,b)
    std::cout << "\n======================================================\n";
    std::cout << "[Part 2] Testing Covariant & Structured Denominators D:\n";
    Poly4 alpha, beta, gamma, delta;
    alpha[Mon4{{0,0,6,0}}] = 7; alpha[Mon4{{0,0,5,1}}] = 5; alpha[Mon4{{0,0,1,5}}] = 15; alpha[Mon4{{0,0,0,6}}] = 7;
    beta[Mon4{{1,1,3,1}}] = 10;
    gamma[Mon4{{3,1,1,1}}] = 5;
    delta[Mon4{{6,0,0,0}}] = 7; delta[Mon4{{5,1,0,0}}] = 5; delta[Mon4{{1,5,0,0}}] = 10; delta[Mon4{{0,6,0,0}}] = 7;

    // Binary cubic Hessian components:
    // h2 = beta^2 - 3*alpha*gamma (bidegree 4,8)
    // h1 = beta*gamma - 9*alpha*delta (bidegree 6,6)
    // h0 = gamma^2 - 3*beta*delta (bidegree 8,4)
    Poly4 h2 = poly_sub(poly_mul(beta, beta), poly_scale(poly_mul(alpha, gamma), 3));
    Poly4 h1 = poly_sub(poly_mul(beta, gamma), poly_scale(poly_mul(alpha, delta), 9));
    Poly4 h0 = poly_sub(poly_mul(gamma, gamma), poly_scale(poly_mul(beta, delta), 3));
    Poly4 H_cub = poly_add(h2, poly_add(h1, h0));

    // Binary cubic Discriminant Delta_cub (degree 24, bidegree 12,12):
    // b^2*g^2 - 4*a*g^3 - 4*b^3*d - 27*a^2*d^2 + 18*a*b*g*d
    Poly4 disc_cub = poly_mul(poly_mul(beta, beta), poly_mul(gamma, gamma));
    disc_cub = poly_sub(disc_cub, poly_scale(poly_mul(alpha, poly_pow(gamma, 3)), 4));
    disc_cub = poly_sub(disc_cub, poly_scale(poly_mul(poly_pow(beta, 3), delta), 4));
    disc_cub = poly_sub(disc_cub, poly_scale(poly_mul(poly_mul(alpha, alpha), poly_mul(delta, delta)), 27));
    disc_cub = poly_add(disc_cub, poly_scale(poly_mul(poly_mul(alpha, beta), poly_mul(gamma, delta)), 18));

    Poly4 pa_var{{Mon4{{1,0,0,0}}, 1}}, pb_var{{Mon4{{0,1,0,0}}, 1}};
    Poly4 pc_var{{Mon4{{0,0,1,0}}, 1}}, pd_var{{Mon4{{0,0,0,1}}, 1}};
    Poly4 a2_2c2 = poly_add(poly_pow(pa_var, 2), poly_scale(poly_pow(pc_var, 2), 2));
    Poly4 abcd = poly_mul(poly_mul(pa_var, pb_var), poly_mul(pc_var, pd_var));

    // 4x4 Hessian determinant of H7 (degree 4*4 = 16)!
    Poly4 Hmat[4][4];
    for (int i=0;i<4;++i) for (int j=0;j<4;++j) Hmat[i][j] = poly_deriv(poly_deriv(H7, i), j);
    auto det3 = [&](int r0, int r1, int r2, int c0, int c1, int c2) {
        int rs[3] = {r0,r1,r2}, cs[3] = {c0,c1,c2};
        Poly4 t0 = poly_mul(Hmat[rs[0]][cs[0]], poly_sub(poly_mul(Hmat[rs[1]][cs[1]], Hmat[rs[2]][cs[2]]), poly_mul(Hmat[rs[1]][cs[2]], Hmat[rs[2]][cs[1]])));
        Poly4 t1 = poly_mul(Hmat[rs[0]][cs[1]], poly_sub(poly_mul(Hmat[rs[1]][cs[0]], Hmat[rs[2]][cs[2]]), poly_mul(Hmat[rs[1]][cs[2]], Hmat[rs[2]][cs[0]])));
        Poly4 t2 = poly_mul(Hmat[rs[0]][cs[2]], poly_sub(poly_mul(Hmat[rs[1]][cs[0]], Hmat[rs[2]][cs[1]]), poly_mul(Hmat[rs[1]][cs[1]], Hmat[rs[2]][cs[0]])));
        return poly_add(poly_sub(t0, t1), t2);
    };
    Poly4 Hess4 = poly_sub(
        poly_add(poly_mul(Hmat[0][0], det3(1,2,3, 1,2,3)), poly_mul(Hmat[0][2], det3(1,2,3, 0,1,3))),
        poly_add(poly_mul(Hmat[0][1], det3(1,2,3, 0,2,3)), poly_mul(Hmat[0][3], det3(1,2,3, 0,1,2)))
    );

    std::vector<std::pair<std::string, Poly4>> cands = {
        {"H7", H7},
        {"H7^2", poly_pow(H7, 2)},
        {"f_a", delta},
        {"f_c", alpha},
        {"f_a * f_c", poly_mul(delta, alpha)},
        {"H7 * f_a * f_c", poly_mul(H7, poly_mul(delta, alpha))},
        {"H_cub (deg 12)", H_cub},
        {"H7 * H_cub (deg 18)", poly_mul(H7, H_cub)},
        {"H7 * a * c * H_cub (deg 20)", poly_mul(poly_mul(H7, poly_mul(pa_var, pc_var)), H_cub)},
        {"h1 = beta*gamma - 9*alpha*delta (deg 12)", h1},
        {"H7 * h1 (deg 18)", poly_mul(H7, h1)},
        {"Hess4(H7) (deg 16)", Hess4},
        {"H7 * Hess4(H7) (deg 22)", poly_mul(H7, Hess4)},
        {"disc_cub (deg 24)", disc_cub},
        {"H7 * disc_cub (deg 30)", poly_mul(H7, disc_cub)},
        {"H7 * a^2 * c^2 * disc_cub (deg 34)", poly_mul(poly_mul(H7, poly_mul(poly_pow(pa_var,2), poly_pow(pc_var,2))), disc_cub)},
        {"H7 * abcd * (a^2+2c^2) (deg 12)", poly_mul(H7, poly_mul(abcd, a2_2c2))},
        {"H_cub * f_a * f_c (deg 24)", poly_mul(H_cub, poly_mul(delta, alpha))},
        {"H7 * H_cub * f_a * f_c (deg 30)", poly_mul(H7, poly_mul(H_cub, poly_mul(delta, alpha)))}
    };
    for (auto& [lbl, pol] : cands) {
        test_fixed_D(lbl, pol, H7, A7);
    }
    return 0;
}
