#include <iostream>
#include <vector>
#include <map>
#include <numeric>
#include <algorithm>
#include <cstdint>

// Reuse audit_WD by including audit_source_image_all.cpp with main renamed
#define main old_main
#include "audit_source_image_all.cpp"
#undef main

static bool is_prime(u64 n) {
    if (n < 2) return false;
    for (u64 d = 2; d * d <= n; ++d) if (n % d == 0) return false;
    return true;
}

int main() {
    Poly4 A_src{{Mon4{{3,0,3,0}}, 1}};
    for (u64 p = 11; p < 400; ++p) {
        if (!is_prime(p) || p == 7) continue;
        ModCtx Mp{p};
        u64 kp = Mp.mul(5, Mp.inv(7));
        auto find_CB = [&](u64 k5, u64 k1) -> std::pair<std::vector<u64>, std::vector<u64>> {
            for (u64 c0 = 0; c0 < p; ++c0)
            for (u64 c1 = 0; c1 < p; ++c1)
            for (u64 c2 = 0; c2 < p; ++c2) {
                u64 C[3] = {c0, c1, c2};
                u64 C2[5] = {0}, C3[7] = {0};
                for (int i=0;i<3;++i) for (int j=0;j<3;++j) C2[i+j] = Mp.add(C2[i+j], Mp.mul(C[i], C[j]));
                for (int i=0;i<5;++i) for (int j=0;j<3;++j) C3[i+j] = Mp.add(C3[i+j], Mp.mul(C2[i], C[j]));
                u64 G[7] = {
                    Mp.sub(C3[0], 1), Mp.sub(C3[1], k5), C3[2], C3[3], C3[4], Mp.sub(C3[5], k1), Mp.sub(C3[6], 1)
                };
                for (u64 b0 = 1; b0 < p; ++b0) {
                    if (Mp.mul(b0, b0) != G[0]) continue;
                    u64 inv2b0 = Mp.inv(Mp.mul(2, b0));
                    u64 b1 = Mp.mul(G[1], inv2b0);
                    u64 b2 = Mp.mul(Mp.sub(G[2], Mp.mul(b1, b1)), inv2b0);
                    u64 b3 = Mp.mul(Mp.sub(G[3], Mp.mul(2, Mp.mul(b1, b2))), inv2b0);
                    if (Mp.add(Mp.mul(b2, b2), Mp.mul(2, Mp.mul(b1, b3))) != G[4]) continue;
                    if (Mp.mul(2, Mp.mul(b2, b3)) != G[5]) continue;
                    if (Mp.mul(b3, b3) != G[6]) continue;
                    return {{c0, c1, c2}, {b0, b1, b2, b3}};
                }
            }
            return {{}, {}};
        };
        auto [Ca, Ba] = find_CB(kp, Mp.mul(2, kp));
        if (Ba.empty()) continue;
        auto [Cc, Bc] = find_CB(kp, Mp.mul(3, kp));
        if (Bc.empty()) continue;

        std::cout << "FOUND SIMULTANEOUS BINARY CONTACTS AT p = " << p << "!\n";
        Poly4 pBa, pBc, pCa, pCc;
        for (int j=0;j<4;++j) { pBa[Mon4{{3-j,j,0,0}}] = Ba[j]; pBc[Mon4{{0,0,3-j,j}}] = Bc[j]; }
        for (int j=0;j<3;++j) { pCa[Mon4{{2-j,j,0,0}}] = Ca[j]; pCc[Mon4{{0,0,2-j,j}}] = Cc[j]; }
        Poly4 D_bin = poly_sub(poly_mul(pBa, pBa, Mp), poly_mul(pBc, pBc, Mp), Mp);
        Poly4 va{{Mon4{{1,0,0,0}}, 1}}, vc{{Mon4{{0,0,1,0}}, 1}};
        Poly4 H0 = make_H(kp, Mp, false);
        Poly4 H1 = make_H(kp, Mp, true);
        int e0 = compute_poly_exact_dim(H0, Mp);
        int e1 = compute_poly_exact_dim(H1, Mp);

        for (auto& [lbl, D] : std::vector<std::pair<std::string, Poly4>>{
            {"D_bin = B_a^2 - B_c^2 (deg 6)", D_bin},
            {"a * D_bin (deg 7)", poly_mul(va, D_bin, Mp)},
            {"a^2 * D_bin (deg 8)", poly_mul(poly_pow(va,2,Mp), D_bin, Mp)},
            {"a*c * D_bin (deg 8)", poly_mul(poly_mul(va,vc,Mp), D_bin, Mp)},
            {"B_a^2 + B_c^2 (deg 6)", poly_add(poly_mul(pBa,pBa,Mp), poly_mul(pBc,pBc,Mp), Mp)},
            {"a^2 * (B_a^2 + B_c^2) (deg 8)", poly_mul(poly_pow(va,2,Mp), poly_add(poly_mul(pBa,pBa,Mp), poly_mul(pBc,pBc,Mp), Mp), Mp)}
        }) {
            auto r0 = audit_WD(D, H0, A_src, Mp, e0);
            auto r1 = audit_WD(D, H1, A_src, Mp, e1);
            std::cout << "  [" << lbl << "] H0: dimW=" << r0.dim_WD << " (new=" << r0.new_classes << ")"
                      << " | H_coupled: dimW=" << r1.dim_WD << " (new=" << r1.new_classes << ")\n";
        }
        break;
    }
    return 0;
}
