/* Exhaustive F_11-rational split-quadric fiber. Compile with c++ -O3
 * -std=c++17. */
#include <array>
#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>
using namespace std;
const int p = 11, k = 7; // 5/7 in F_11
int reduce(int a) {
  int b = a % p;
  return b < 0 ? b + p : b;
}
int power(int a, int n) {
  int v = 1;
  while (n) {
    if (n % 2)
      v = v * a % p;
    a = a * a % p;
    n /= 2;
  }
  return v;
}
int inverse(int a) {
  a = reduce(a);
  assert(a);
  return power(a, 9);
}
int fx(int a) { return reduce(1 + power(a, 6) + k * (a + 2 * power(a, 5))); }
int fy(int a) { return reduce(1 + power(a, 6) + k * (a + 3 * power(a, 5))); }
using Cubic = array<int, 4>;
using Poly = array<array<int, 7>, 7>;
int index4(int a, int b, int c, int d) { return ((a * p + b) * p + c) * p + d; }
vector<Cubic> axis[p * p * p * p];
int edges[p][p][p][p];
void init() {
  for (int a = 1; a <= p / 2; a++)
    for (int b = 0; b < p; b++)
      for (int c = 0; c < p; c++)
        for (int d = 0; d < p; d++)
          axis[index4(a * a % p, reduce(2 * a * c - b * b),
                      reduce(c * c - 2 * b * d), reduce(-d * d))]
              .push_back({a, b, c, d});
  for (int i = 0; i < p; i++)
    for (int j = 0; j < p; j++)
      for (int A = 0; A < p; A++)
        for (int B = 0; B < p; B++)
          edges[i][j][A][B] = index4(
              reduce(fx(i) - A * A * A), reduce(k * i * j - 3 * A * A * B),
              reduce(2 * k * i * j - 3 * A * B * B), reduce(fy(j) - B * B * B));
}
Poly multiply(const Poly &a, const Poly &b, int n, int m) {
  Poly z{};
  for (int i = 0; i <= n; i++)
    for (int j = 0; j <= n; j++)
      if (a[i][j])
        for (int h = 0; h <= m; h++)
          for (int l = 0; l <= m; l++)
            z[i + h][j + l] += a[i][j] * b[h][l];
  for (int i = 0; i <= n + m; i++)
    for (int j = 0; j <= n + m; j++)
      z[i][j] %= p;
  return z;
}
void add(Poly &z, const Poly &a, int c, int n) {
  for (int i = 0; i <= n; i++)
    for (int j = 0; j <= n; j++)
      z[i][j] = reduce(z[i][j] + a[i][j] * c);
}
void chart(int m, int n, int u, int v, int K, Poly &H, array<Poly, 5> &base) {
  array<array<Poly, 7>, 4> f{};
  f[0][1][0][0] = 1;
  f[0][1][1][1] = K;
  f[1][1][0][0] = n;
  f[1][1][1][1] = m * K % p;
  f[2][1][1][0] = f[2][1][0][1] = 1;
  f[3][1][1][0] = v;
  f[3][1][0][1] = u;
  for (int i = 0; i < 4; i++)
    for (int j = 2; j <= 6; j++)
      f[i][j] = multiply(f[i][j - 1], f[i][1], j - 1, 1);
  for (int i = 0; i < 4; i++)
    add(H, f[i][6], 1, 6);
  add(H, multiply(f[0][5], f[1][1], 5, 1), k, 6);
  add(H, multiply(f[1][5], f[0][1], 5, 1), 2 * k, 6);
  add(H, multiply(f[2][5], f[3][1], 5, 1), k, 6);
  add(H, multiply(f[3][5], f[2][1], 5, 1), 3 * k, 6);
  auto ab = multiply(f[0][1], f[1][1], 1, 1);
  auto cd = multiply(f[2][1], f[3][1], 1, 1);
  auto abcd = multiply(ab, cd, 2, 2);
  add(H, multiply(abcd, f[0][2], 4, 2), k, 6);
  add(H, multiply(abcd, f[2][2], 4, 2), 2 * k, 6);
  base = {f[0][2], f[1][2], f[2][2], cd, f[3][2]};
}
int normcoeff(const Poly &h, int a, int b) {
  int t = 0;
  for (int i = max(0, a - 3); i <= min(3, a); i++)
    for (int j = max(0, b - 3); j <= min(3, b); j++)
      t += ((a + b - i - j) % 2 ? -1 : 1) * h[i][j] * h[a - i][b - j];
  return reduce(t);
}
bool finish(const Poly &f, const Cubic &z, const Cubic &w, Poly &h) {
  auto [u, v, y, x] = z;
  if (!x || !reduce(u * x - v * y))
    return false;
  int iv = inverse(u * x - v * y);
  for (int j = 0; j < 4; j++)
    h[0][j] = z[j];
  for (int i = 1; i <= 3; i++) {
    auto half = [&](int j) {
      return reduce((f[i][j] - normcoeff(h, i, j)) * 6);
    };
    int c0 = w[i], c1, c2, c3;
    if (i % 2) {
      c1 = reduce((half(1) + v * c0) * inverse(u));
      int b = reduce(half(3) - y * c1 + x * c0), d = half(5);
      c2 = reduce((y * b - u * d) * iv);
      c3 = reduce((x * b - v * d) * iv);
    } else {
      c3 = reduce(-half(6) * inverse(x));
      int b = reduce(half(2) - y * c0), d = reduce(half(4) + v * c3);
      c1 = reduce((y * b - u * d) * iv);
      c2 = reduce((x * b - v * d) * iv);
    }
    h[i][0] = c0;
    h[i][1] = c1;
    h[i][2] = c2;
    h[i][3] = c3;
  }
  for (int i = 4; i <= 6; i++)
    for (int j = 1 + (i % 2 == 0); j <= 6; j += 2)
      if (normcoeff(h, i, j) != f[i][j])
        return false;
  for (int i = 0; i <= 6; i++)
    for (int j = 0; j <= 6; j++)
      assert(normcoeff(h, i, j) == f[i][j]);
  return true;
}
int main(int argc, char **argv) {
  bool census = argc > 1 && string(argv[1]) == "--census";
  int r = argc > 1 ? reduce(atoi(argv[1])) : 4,
      s = argc > 2 ? reduce(atoi(argv[2])) : 6;
  init();
  long long quads = 0, corners = 0, axis_pairs = 0, contacts = 0, points = 0;
  for (int m = 0; m < p; m++)
    for (int n = m + 1; n < p; n++) {
      if (!reduce(m + n) || (!census && (r == m || r == n)))
        continue;
      int D = reduce(-inverse(m + n));
      for (int u = 0; u < p; u++)
        for (int v = u + 1; v < p; v++) {
          if (census) {
            s = 0;
            while (s == u || s == v)
              s++;
          } else if (s == u || s == v)
            continue;
          vector<array<int, 6>> seeds;
          for (int M = 0; M < p; M++)
            for (int N = 0; N < p; N++) {
              int Q = reduce((N - M) * inverse(n * n - m * m)),
                  P = reduce(M - Q * m * m);
              for (int Y = 0; Y < p; Y++)
                if (!axis[edges[m][u][M][Y]].empty() &&
                    !axis[edges[n][u][N][Y]].empty())
                  for (int Z = 0; Z < p; Z++)
                    if (!axis[edges[m][v][M][Z]].empty() &&
                        !axis[edges[n][v][N][Z]].empty())
                      seeds.push_back({P, Q, N, Y, Z, M});
            }
          for (int e = 1; e < p; e++) {
            int t = census ? 1
                           : reduce(-D * (r - m) * (r - n) *
                                    inverse(e * (s - u) * (s - v)));
            if (power(t, 5) != 1)
              continue;
            quads++;
            corners += seeds.size();
            int Hpt = reduce(fx(r) + k * r * s * t + 2 * k * r * s * t * t +
                             fy(s) * t * t * t);
            int l = power(Hpt, 7), K = reduce(-e * (v - u) * (v - u) *
                                              inverse(D * (n - m) * (n - m)));
            assert(power(l, 3) == Hpt);
            Poly H{};
            array<Poly, 5> base;
            chart(m, n, u, v, K, H, base);
            for (auto seed : seeds) {
              auto [P, Q, N, Y, Z, M] = seed;
              for (int lift = 0; lift < (census ? p : 1); lift++) {
                int W =
                    census ? lift : reduce((l - P - r * r * Q) * inverse(t));
                auto co = [&](int a, int b, int c, int y) {
                  return reduce(y * inverse((c - a) * (c - b)));
                };
                int iu = co(v, s, u, Y), iv = co(u, s, v, Z),
                    is = co(u, v, s, W);
                int R = reduce(iu * v * s + iv * u * s + is * u * v),
                    T = reduce(-iu * (v + s) - iv * (u + s) - is * (u + v)),
                    S = reduce(iu + iv + is);
                array<int, 5> vals = {P, Q, R, T, S};
                Poly q{};
                for (int j = 0; j < 5; j++)
                  add(q, base[j], vals[j], 2);
                Poly f = H;
                add(f, multiply(multiply(q, q, 2, 2), q, 4, 2), -1, 6);
                for (auto z : axis[edges[n][u][N][Y]])
                  for (auto w : axis[edges[n][v][N][Z]]) {
                    if (!z[3] || !reduce(z[0] * z[3] - z[1] * z[2]))
                      continue;
                    axis_pairs++;
                    Poly h{};
                    if (!finish(f, z, w, h))
                      continue;
                    contacts++;
                    if (census) {
                      cout << "contact " << m << " " << n << " " << u << " "
                           << v << " " << e << " Q1 " << P << " " << Q << " "
                           << R << " " << T << " " << S << " h";
                      for (int i = 0; i < 4; i++)
                        for (int j = 0; j < 4; j++)
                          cout << " " << h[i][j];
                      cout << endl;
                      continue;
                    }
                    for (int a = 1; a < p; a++)
                      if (a * a % p == t) {
                        int A = reduce((m - n) * inverse(m - r)),
                            ww = reduce(a * A * (s - u) * inverse(v - u)),
                            zz = reduce(a * A * (v - s) * inverse(v - u));
                        assert(reduce(1 + K * ww * zz) == A &&
                               reduce(n + m * K * ww * zz) == r * A % p);
                        assert(reduce(ww + zz) == a * A % p &&
                               reduce(v * ww + u * zz) == s * a * A % p);
                        int hp = 0;
                        for (int i = 0; i < 4; i++)
                          for (int j = 0; j < 4; j++)
                            hp += h[i][j] * power(ww, i) * power(zz, j);
                        if (reduce(hp))
                          continue;
                        points++;
                        cout << "point roots " << m << " " << n << " " << u
                             << " " << v << " scale " << e << " Q1 " << P << " "
                             << Q << " " << R << " " << T << " " << S
                             << " alpha " << a << " h";
                        for (int i = 0; i < 4; i++)
                          for (int j = 0; j < 4; j++)
                            cout << " " << h[i][j];
                        cout << "\n";
                      }
                  }
              }
            }
          }
        }
    }
  cout << "ratios " << r << " " << s << " quadrics " << quads
       << " linear-candidates " << quads * p * p * p * p * (census ? p : 1)
       << " four-corners " << corners << " axis-pairs " << axis_pairs
       << " contacts-mod-sign " << contacts << " points-mod-sign " << points
       << "\n";
}
