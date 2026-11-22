// Scalar-discrepancy approximant basis over F_11. Input m rows of n bytes,
// row j = coefficients of Z(x)^j.  Optional affine shift s_j=j*weight.
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
#if defined(__ARM_NEON)
#include <arm_neon.h>
#endif
using Byte = std::uint8_t;
Byte submul[11][16];
void subtract(Byte *row, const Byte *pivot, int n, int coeff) {
  int j = 0;
#if defined(__ARM_NEON)
  auto tab = vld1q_u8(submul[coeff]);
  auto prime = vdupq_n_u8(11);
  for (; j + 16 <= n; j += 16) {
    auto z = vld1q_u8(pivot + j);
    auto a = vld1q_u8(row + j);
    auto s = vsubq_u8(vaddq_u8(a, prime), vqtbl1q_u8(tab, z));
    vst1q_u8(row + j, vsubq_u8(s, vandq_u8(vcgeq_u8(s, prime), prime)));
  }
#endif
  for (; j < n; j++) {
    int z = row[j] + 11 - submul[coeff][pivot[j]];
    row[j] = z >= 11 ? z - 11 : z;
  }
}
int main(int argc, char **argv) {
  if (argc != 5) {
    std::cerr << "usage: sweep_pade powers.bin m n weight\n";
    return 2;
  }
  const int m = std::atoi(argv[2]), n = std::atoi(argv[3]),
            weight = std::atoi(argv[4]);
  if (m < 2 || n < 1 || weight < 0)
    return 2;
  for (int i = 0; i < 11; i++)
    for (int j = 0; j < 11; j++)
      submul[i][j] = i * j % 11;
  int inverse[11] = {0, 1, 6, 4, 3, 9, 2, 8, 7, 5, 10};
  std::vector<Byte> a(static_cast<size_t>(m) * n);
  std::ifstream in(argv[1], std::ios::binary);
  in.read(reinterpret_cast<char *>(a.data()), a.size());
  if (!in)
    return 3;
  auto row = [&](int i) { return a.data() + static_cast<size_t>(i) * n; };
  std::vector<int> degree(m);
  for (int i = 0; i < m; i++)
    degree[i] = weight * i;
  for (int k = 0; k < n; k++) {
    int pivot = -1;
    for (int i = 0; i < m; i++)
      if (row(i)[k] && (pivot < 0 || degree[i] < degree[pivot]))
        pivot = i;
    if (pivot >= 0) {
      Byte *p = row(pivot);
      int inv = inverse[p[k]];
      for (int i = 0; i < m; i++)
        if (i != pivot && row(i)[k])
          subtract(row(i) + k, p + k, n - k, row(i)[k] * inv % 11);
      std::memmove(p + k + 1, p + k, n - k - 1);
      p[k] = 0;
      degree[pivot]++;
    }
  }
  for (int i = 0; i < m; i++)
    std::cout << (i ? "," : "") << degree[i];
  std::cout << "\n";
}
