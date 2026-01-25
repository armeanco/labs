//https://judge.yosupo.jp/problem/many_factorials
#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <tuple>
#include <stack>
#include <queue>
#include <deque>
#include <algorithm>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <bitset>
#include <cmath>
#include <functional>
#include <cassert>
#include <climits>
#include <iomanip>
#include <numeric>
#include <memory>
#include <random>
#include <thread>
#include <chrono>
#define allof(obj) (obj).begin(), (obj).end()
#define range(i, l, r) for(int i=l;i<r;i++)
#define unique_elem(obj) obj.erase(std::unique(allof(obj)), obj.end())
#define bit_subset(i, S) for(int i=S, zero_cnt=0;(zero_cnt+=i==S)<2;i=(i-1)&S)
#define bit_kpop(i, n, k) for(int i=(1<<k)-1,x_bit,y_bit;i<(1<<n);x_bit=(i&-i),y_bit=i+x_bit,i=(!i?(1<<n):((i&~y_bit)/x_bit>>1)|y_bit))
#define bit_kth(i, k) ((i >> k)&1)
#define bit_highest(i) (i?63-__builtin_clzll(i):-1)
#define bit_lowest(i) (i?__builtin_ctzll(i):-1)
#define sleepms(t) std::this_thread::sleep_for(std::chrono::milliseconds(t))
using ll = long long;
using ld = long double;
using ul = uint64_t;
using pi = std::pair<int, int>;
using pl = std::pair<ll, ll>;
using namespace std;

template<typename F, typename S>
std::ostream &operator << (std::ostream &dest, const std::pair<F, S> &p) {
    dest << p.first << ' ' << p.second;
    return dest;
}

template<typename A, typename B>
std::ostream &operator << (std::ostream &dest, const std::tuple<A, B> &t) {
    dest << std::get<0>(t) << ' ' << std::get<1>(t);
    return dest;
}

template<typename A, typename B, typename C>
std::ostream &operator << (std::ostream &dest, const std::tuple<A, B, C> &t) {
    dest << std::get<0>(t) << ' ' << std::get<1>(t) << ' ' << std::get<2>(t);
    return dest;
}

template<typename A, typename B, typename C, typename D>
std::ostream &operator << (std::ostream &dest, const std::tuple<A, B, C, D> &t) {
    dest << std::get<0>(t) << ' ' << std::get<1>(t) << ' ' << std::get<2>(t) << ' ' << std::get<3>(t);
    return dest;
}

template<typename T>
std::ostream &operator << (std::ostream &dest, const std::vector<std::vector<T>> &v) {
    int sz = v.size();
    if (!sz) return dest;
    for (int i = 0; i < sz; i++) {
        int m = v[i].size();
        for (int j = 0; j < m; j++) dest << v[i][j] << (i != sz - 1 && j == m - 1 ? '\n' : ' ');
    }
    return dest;
}

template<typename T>
std::ostream &operator << (std::ostream &dest, const std::vector<T> &v) {
    int sz = v.size();
    if (!sz) return dest;
    for (int i = 0; i < sz - 1; i++) dest << v[i] << ' ';
    dest << v[sz - 1];
    return dest;
}

template<typename T, size_t sz>
std::ostream &operator << (std::ostream &dest, const std::array<T, sz> &v) {
    if (!sz) return dest;
    for (int i = 0; i < sz - 1; i++) dest << v[i] << ' ';
    dest << v[sz - 1];
    return dest;
}

template<typename T>
std::ostream &operator << (std::ostream &dest, const std::set<T> &v) {
    for (auto itr = v.begin(); itr != v.end();) {
        dest << *itr;
        itr++;
        if (itr != v.end()) dest << ' ';
    }
    return dest;
}

template<typename T, typename E>
std::ostream &operator << (std::ostream &dest, const std::map<T, E> &v) {
    for (auto itr = v.begin(); itr != v.end(); ) {
        dest << '(' << itr->first << ", " << itr->second << ')';
        itr++;
        if (itr != v.end()) dest << '\n';
    }
    return dest;
}

template<typename T>
vector<T> make_vec(size_t sz, T val) { return std::vector<T>(sz, val); }

template<typename T, typename... Tail>
auto make_vec(size_t sz, Tail ...tail) {
    return std::vector<decltype(make_vec<T>(tail...))>(sz, make_vec<T>(tail...));
}

template<typename T>
vector<T> read_vec(size_t sz) {
    std::vector<T> v(sz);
    for (int i = 0; i < (int)sz; i++) std::cin >> v[i];
    return v;
}

template<typename T, typename... Tail>
auto read_vec(size_t sz, Tail ...tail) {
    auto v = std::vector<decltype(read_vec<T>(tail...))>(sz);
    for (int i = 0; i < (int)sz; i++) v[i] = read_vec<T>(tail...);
    return v;
}

// x / yä»¥ä¸Šã®æœ€å°ã®æ•´æ•°
ll ceil_div(ll x, ll y) {
    assert(y > 0);
    return (x + (x > 0 ? y - 1 : 0)) / y;
}

// x / yä»¥ä¸‹ã®æœ€å¤§ã®æ•´æ•°
ll floor_div(ll x, ll y) {
    assert(y > 0);
    return (x + (x > 0 ? 0 : -y + 1)) / y;
}

void io_init() {
    std::cin.tie(nullptr);
    std::ios::sync_with_stdio(false);
}

#include <type_traits>
#include <fstream>

// @param m `1 <= m`
constexpr long long safe_mod(long long x, long long m){
  x %= m;
  if (x < 0) x += m;
  return x;
}

// x^n mod m
// @param n `0 <= n`
// @param m `1 <= m`
constexpr long long pow_mod_constexpr(long long x, long long n, int m) {
    if (m == 1) return 0;
    unsigned int _m = (unsigned int)(m);
    unsigned long long r = 1;
    unsigned long long y = safe_mod(x, m);
    while (n) {
        if (n & 1) r = (r * y) % _m;
        y = (y * y) % _m;
        n >>= 1;
    }
    return r;
}

constexpr __uint128_t pow_mod64_constexpr(__int128_t x, __uint128_t n, unsigned long long m) {
    if (m == 1) return 0;
    __uint128_t r = 1;
    if (x >= m) x %= m;
    if (x < 0) x += m;
    while (n) {
        if (n & 1) r = (r * x) % m;
        x = (x * x) % m;
        n >>= 1;
    }
    return r;
}

constexpr bool miller_rabin32_constexpr(int n) {
    if (n <= 1) return false;
    if (n == 2 || n == 7 || n == 61) return true;
    if (n % 2 == 0) return false;
    long long d = n - 1;
    while (d % 2 == 0) d /= 2;
    constexpr long long bases[3] = {2, 7, 61};
    for (long long a : bases) {
        long long t = d;
        long long y = pow_mod_constexpr(a, t, n);
        while (t != n - 1 && y != 1 && y != n - 1) { 
            y = y * y % n;
            t <<= 1;
        }
        if (y != n - 1 && t % 2 == 0) {
            return false;
        }
    }
    return true;
}

template<int n>
constexpr bool miller_rabin32 = miller_rabin32_constexpr(n);

// -10^18 <= _a, _b <= 10^18
long long gcd(long long _a, long long _b) {
    long long a = abs(_a), b = abs(_b);
    if (a == 0) return b;
    if (b == 0) return a;
    int shift = __builtin_ctzll(a | b);
    a >>= __builtin_ctzll(a);
    do{
        b >>= __builtin_ctzll(b);
        if(a > b) std::swap(a, b);
        b -= a;
    } while (b);
    return a << shift;
}

// æœ€å¤§ã§a*b
// -10^18 <= a, b <= 10^18
// a, bã¯è² ã§ã‚‚ã„ã„ãŒéžè² ã®å€¤ã‚’è¿”ã™
__int128_t lcm(long long a, long long b) {
    a = abs(a), b = abs(b);
    long long g = gcd(a, b);
    if (!g) return 0;
    return __int128_t(a) * b / g;
}

// {x, y, gcd(a, b)} s.t. ax + by = gcd(a, b)
// g >= 0
std::tuple<long long, long long, long long> extgcd(long long a, long long b) {
    long long x, y;
    for (long long u = y = 1, v = x = 0; a;) {
        long long q = b / a;
        std::swap(x -= q * u, u);
        std::swap(y -= q * v, v);
        std::swap(b -= q * a, a);
    }
    // x + k * (b / g), y - k * (a / g) ã‚‚æ¡ä»¶ã‚’æº€ãŸã™(kã¯ä»»æ„ã®æ•´æ•°)
    return {x, y, b};
}

// @param b `1 <= b`
// @return pair(g, x) s.t. g = gcd(a, b), xa = g (mod b), 0 <= x < b/g
constexpr std::pair<long long, long long> inv_gcd(long long a, long long b) {
    a = safe_mod(a, b);
    if (a == 0) return {b, 0};
    long long s = b, t = a;
    long long m0 = 0, m1 = 1;
    while (t) {
        long long u = s / t;
        s -= t * u;
        m0 -= m1 * u;
        auto tmp = s;
        s = t;
        t = tmp;
        tmp = m0;
        m0 = m1;
        m1 = tmp;
    }
    if (m0 < 0) m0 += b / s;
    return {s, m0};
}

template <int m, std::enable_if_t<(1 <= m)>* = nullptr>
struct modint32_static {
    using mint = modint32_static;
  public:
    static constexpr int mod() { return m; }
    
    static mint raw(int v) {
        mint x;
        x._v = v;
        return x;
    }
  
    modint32_static(): _v(0) {}
    
    template <class T>
    modint32_static(T v) { 
        long long x = v % (long long)umod();
        if (x < 0) x += umod();
        _v = x;
    }

    unsigned int val() const { return _v; }
    
    mint& operator ++ () {
        _v++;
        if (_v == umod()) _v = 0;
        return *this;
    }
    mint& operator -- () {
        if (_v == 0) _v = umod();
        _v--;
        return *this;
    }
    mint operator ++ (int) {
        mint result = *this;
        ++*this;
        return result;
    }
    mint operator -- (int) {
        mint result = *this;
        --*this;
        return result;
    }
    mint& operator += (const mint& rhs) {
        _v += rhs._v;
        if (_v >= umod()) _v -= umod();
        return *this;
    }
    mint& operator -= (const mint& rhs) {
        _v -= rhs._v;
        if (_v >= umod()) _v += umod();
        return *this;
    }
    mint& operator *= (const mint& rhs) {
        unsigned long long z = _v;
        z *= rhs._v;
        _v = (unsigned int)(z % umod());
        return *this;
    }
    mint& operator /= (const mint& rhs) { return *this = *this * rhs.inv(); }
    mint operator + () const { return *this; }
    mint operator - () const { return mint() - *this; }
    mint pow(long long n) const {
        assert(0 <= n);
        mint x = *this, r = 1;
        while (n) {
            if (n & 1) r *= x;
            x *= x;
            n >>= 1;
        }
        return r;
    }
    mint inv() const {
        if (prime) {
            assert(_v);
            return pow(umod() - 2);
        } else {
            auto eg = inv_gcd(_v, m);
            assert(eg.first == 1);
            return eg.second;
        }
    }
    friend mint operator + (const mint& lhs, const mint& rhs) { return mint(lhs) += rhs; }
    friend mint operator - (const mint& lhs, const mint& rhs) { return mint(lhs) -= rhs; }
    friend mint operator * (const mint& lhs, const mint& rhs) { return mint(lhs) *= rhs; }
    friend mint operator / (const mint& lhs, const mint& rhs) { return mint(lhs) /= rhs; }
    friend bool operator == (const mint& lhs, const mint& rhs) { return lhs._v == rhs._v; }
    friend bool operator != (const mint& lhs, const mint& rhs) { return lhs._v != rhs._v; }
  private:
    unsigned int _v;
    static constexpr unsigned int umod() { return m; }
    static constexpr bool prime = miller_rabin32<m>;
};

template<int m>
std::ostream &operator<<(std::ostream &dest, const modint32_static<m> &a) {
    dest << a.val();
    return dest;
}

using modint998244353 = modint32_static<998244353>;
using modint1000000007 = modint32_static<1000000007>;

// {x^2 â‰¡ aã¨ãªã‚‹xãŒå­˜åœ¨ã™ã‚‹ã‹, x}
// ç´ æ•°mod
template<typename mint>
std::pair<bool, mint> sqrt_mod(mint a) {
    static std::random_device seed_gen;
    static std::mt19937_64 engine(seed_gen());

	if (a == 0) return {true, 0};
	if (mint::mod() == 2) return {true, 1};
    
    // ã‚ªã‚¤ãƒ©ãƒ¼ã®è¦æº–
    if (a.pow((mint::mod() - 1) / 2) != 1) {
        return {false, 0};
    }
    
	if (mint::mod() % 4 == 3) {
        return {true, a.pow(mint::mod() / 4 + 1)};
    }

	long long q = mint::mod() - 1, m = 0;
	while (q % 2 == 0) q >>= 1, m++;
	
	mint z;
    while (true) {
       z = engine();
       if (z.pow((mint::mod() - 1) / 2) == -1) break;
	}
	
  mint c = z.pow(q);
	mint t = a.pow(q);
	mint r = a.pow((q + 1) / 2);
	
    for(; m > 1; m--) {
        if (t.pow(1LL << (m - 2)) != 1) {
            r *= c;
            t *= c * c;
        }
		c *= c;
	}
	return {true, r};
}

constexpr int primitive_root32_constexpr(int m) {
    if (m == 2) return 1;
    if (m == 167772161) return 3;
    if (m == 469762049) return 3;
    if (m == 754974721) return 11;
    if (m == 998244353) return 3;
    int divs[20] = {};
    divs[0] = 2;
    int cnt = 1;
    int x = (m - 1) / 2;
    while (x % 2 == 0) x /= 2;
    for (int i = 3; (long long)(i)*i <= x; i += 2) {
        if (x % i == 0) {
            divs[cnt++] = i;
            while (x % i == 0) {
                x /= i;
            }
        }
    }
    if (x > 1) divs[cnt++] = x;
    for (int g = 2;; g++) {
        bool ok = true;
        for (int i = 0; i < cnt; i++) {
            if (pow_mod_constexpr(g, (m - 1) / divs[i], m) == 1) {
                ok = false;
                break;
            }
        }
        if (ok) return g;
    }
}

template <int m>
constexpr int primitive_root32 = primitive_root32_constexpr(m);

constexpr unsigned int bit_ceil(unsigned int n) {
    unsigned int x = 1;
    while (x < (unsigned int)(n)) x *= 2;
    return x;
}

constexpr int bit_ceil_log(unsigned int n) {
    int x = 0;
    while ((1 << x) < (unsigned int)(n)) x++;
    return x;
}

template <class mint, int g = primitive_root32<mint::mod()>>
struct fft_info {
    static constexpr int rank2 = __builtin_ctz(mint::mod() - 1);
    std::array<mint, rank2 + 1> root;   // root[i]^(2^i) == 1
    std::array<mint, rank2 + 1> iroot;  // root[i] * iroot[i] == 1

    std::array<mint, std::max(0, rank2 - 2 + 1)> rate2;
    std::array<mint, std::max(0, rank2 - 2 + 1)> irate2;

    std::array<mint, std::max(0, rank2 - 3 + 1)> rate3;
    std::array<mint, std::max(0, rank2 - 3 + 1)> irate3;

    fft_info() {
        root[rank2] = mint(g).pow((mint::mod() - 1) >> rank2);
        iroot[rank2] = root[rank2].inv();
        for (int i = rank2 - 1; i >= 0; i--) {
            root[i] = root[i + 1] * root[i + 1];
            iroot[i] = iroot[i + 1] * iroot[i + 1];
        }

        {
            mint prod = 1, iprod = 1;
            for (int i = 0; i <= rank2 - 2; i++) {
                rate2[i] = root[i + 2] * prod;
                irate2[i] = iroot[i + 2] * iprod;
                prod *= iroot[i + 2];
                iprod *= root[i + 2];
            }
        }
        {
            mint prod = 1, iprod = 1;
            for (int i = 0; i <= rank2 - 3; i++) {
                rate3[i] = root[i + 3] * prod;
                irate3[i] = iroot[i + 3] * iprod;
                prod *= iroot[i + 3];
                iprod *= root[i + 3];
            }
        }
    }
};

template <class mint>
void butterfly(std::vector<mint>& a) {
    int n = int(a.size());
    int h = __builtin_ctz((unsigned int)n);

    static const fft_info<mint> info;
    int len = 0;  // a[i, i+(n>>len), i+2*(n>>len), ..] is transformed
    while (len < h) {
        if (h - len == 1) {
            int p = 1 << (h - len - 1);
            mint rot = 1;
            for (int s = 0; s < (1 << len); s++) {
                int offset = s << (h - len);
                for (int i = 0; i < p; i++) {
                    auto l = a[i + offset];
                    auto r = a[i + offset + p] * rot;
                    a[i + offset] = l + r;
                    a[i + offset + p] = l - r;
                }
                if (s + 1 != (1 << len)) {
                    rot *= info.rate2[__builtin_ctz(~(unsigned int)(s))];
                }
            }
            len++;
        } else {
            // 4-base
            int p = 1 << (h - len - 2);
            mint rot = 1, imag = info.root[2];
            for (int s = 0; s < (1 << len); s++) {
                mint rot2 = rot * rot;
                mint rot3 = rot2 * rot;
                int offset = s << (h - len);
                for (int i = 0; i < p; i++) {
                    auto mod2 = 1ULL * mint::mod() * mint::mod();
                    auto a0 = 1ULL * a[i + offset].val();
                    auto a1 = 1ULL * a[i + offset + p].val() * rot.val();
                    auto a2 = 1ULL * a[i + offset + 2 * p].val() * rot2.val();
                    auto a3 = 1ULL * a[i + offset + 3 * p].val() * rot3.val();
                    auto a1na3imag = 1ULL * mint(a1 + mod2 - a3).val() * imag.val();
                    auto na2 = mod2 - a2;
                    a[i + offset] = a0 + a2 + a1 + a3;
                    a[i + offset + 1 * p] = a0 + a2 + (2 * mod2 - (a1 + a3));
                    a[i + offset + 2 * p] = a0 + na2 + a1na3imag;
                    a[i + offset + 3 * p] = a0 + na2 + (mod2 - a1na3imag);
                }
                if (s + 1 != (1 << len)) {
                    rot *= info.rate3[__builtin_ctz(~(unsigned int)(s))];
                }
            }
            len += 2;
        }
    }
}

template <class mint>
void butterfly_inv(std::vector<mint>& a) {
    int n = int(a.size());
    int h = __builtin_ctz((unsigned int)n);

    static const fft_info<mint> info;

    int len = h;  // a[i, i+(n>>len), i+2*(n>>len), ..] is transformed
    while (len) {
        if (len == 1) {
            int p = 1 << (h - len);
            mint irot = 1;
            for (int s = 0; s < (1 << (len - 1)); s++) {
                int offset = s << (h - len + 1);
                for (int i = 0; i < p; i++) {
                    auto l = a[i + offset];
                    auto r = a[i + offset + p];
                    a[i + offset] = l + r;
                    a[i + offset + p] = (unsigned long long)(mint::mod() + l.val() - r.val()) * irot.val();
                }
                if (s + 1 != (1 << (len - 1))) {
                    irot *= info.irate2[__builtin_ctz(~(unsigned int)(s))];
                }
            }
            len--;
        } else {
            // 4-base
            int p = 1 << (h - len);
            mint irot = 1, iimag = info.iroot[2];
            for (int s = 0; s < (1 << (len - 2)); s++) {
                mint irot2 = irot * irot;
                mint irot3 = irot2 * irot;
                int offset = s << (h - len + 2);
                for (int i = 0; i < p; i++) {
                    auto a0 = 1ULL * a[i + offset + 0 * p].val();
                    auto a1 = 1ULL * a[i + offset + 1 * p].val();
                    auto a2 = 1ULL * a[i + offset + 2 * p].val();
                    auto a3 = 1ULL * a[i + offset + 3 * p].val();
                    auto a2na3iimag = 1ULL * mint((mint::mod() + a2 - a3) * iimag.val()).val();

                    a[i + offset] = a0 + a1 + a2 + a3;
                    a[i + offset + 1 * p] = (a0 + (mint::mod() - a1) + a2na3iimag) * irot.val();
                    a[i + offset + 2 * p] = (a0 + a1 + (mint::mod() - a2) + (mint::mod() - a3)) * irot2.val();
                    a[i + offset + 3 * p] = (a0 + (mint::mod() - a1) + (mint::mod() - a2na3iimag)) * irot3.val();
                }
                if (s + 1 != (1 << (len - 2))) {
                    irot *= info.irate3[__builtin_ctz(~(unsigned int)(s))];
                }
            }
            len -= 2;
        }
    }
}

template <class mint>
std::vector<mint> convolution_naive(const std::vector<mint>& a, const std::vector<mint>& b) {
    int n = int(a.size()), m = int(b.size());
    std::vector<mint> ans(n + m - 1);
    if (n < m) {
        for (int j = 0; j < m; j++) {
            for (int i = 0; i < n; i++) {
                ans[i + j] += a[i] * b[j];
            }
        }
    } else {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                ans[i + j] += a[i] * b[j];
            }
        }
    }
    return ans;
}

template <class mint>
std::vector<mint> convolution_fft(std::vector<mint> a, std::vector<mint> b) {
    int n = int(a.size()), m = int(b.size());
    int z = (int)bit_ceil((unsigned int)(n + m - 1));
    a.resize(z);
    butterfly(a);
    b.resize(z);
    butterfly(b);
    for (int i = 0; i < z; i++) {
        a[i] *= b[i];
    }
    butterfly_inv(a);
    a.resize(n + m - 1);
    mint iz = mint(z).inv();
    for (int i = 0; i < n + m - 1; i++) a[i] *= iz;
    return a;
}

template <class mint>
std::vector<mint> convolution_mod(std::vector<mint>&& a, std::vector<mint>&& b) {
    int n = int(a.size()), m = int(b.size());
    if (!n || !m) return {};
    int z = (int)bit_ceil((unsigned int)(n + m - 1));
    assert((mint::mod() - 1) % z == 0);
    if (std::min(n, m) <= 60) return convolution_naive(std::move(a), std::move(b));
    return convolution_fft(std::move(a), std::move(b));
}
template <class mint>
std::vector<mint> convolution_mod(const std::vector<mint>& a, const std::vector<mint>& b) {
    int n = int(a.size()), m = int(b.size());
    if (!n || !m) return {};
    int z = (int)bit_ceil((unsigned int)(n + m - 1));
    assert((mint::mod() - 1) % z == 0);
    if (std::min(n, m) <= 60) return convolution_naive(a, b);
    return convolution_fft(a, b);
}

template <unsigned int mod = 16661, class T>
std::vector<T> convolution_mod(const std::vector<T>& a, const std::vector<T>& b) {
    int n = int(a.size()), m = int(b.size());
    if (!n || !m) return {};
    using mint = modint32_static<mod>;
    int z = (int)bit_ceil((unsigned int)(n + m - 1));
    assert((mint::mod() - 1) % z == 0);

    std::vector<mint> a2(n), b2(m);
    for (int i = 0; i < n; i++) {
        a2[i] = mint(a[i]);
    }
    for (int i = 0; i < m; i++) {
        b2[i] = mint(b[i]);
    }
    auto c2 = convolution_mod(std::move(a2), std::move(b2));
    std::vector<T> c(n + m - 1);
    for (int i = 0; i < n + m - 1; i++) {
        c[i] = c2[i].val();
    }
    return c;
}

template<typename mint>
struct fps_base {
    // inv(2^k)
    static mint ipow2(int k) {
        static constexpr int max_size = 30;
        assert(k < max_size);
        static std::vector<mint> table;
        if (table.empty()) {
            mint i2 = mint(2).inv();
            table.resize(max_size, 1);
            for (int i = 1; i < max_size; i++) {
                table[i] = table[i - 1] * i2;
            }
        }
        return table[k];
    }

    static mint inv_low(int i) {
        static std::vector<mint> table{0};
        if (table.size() <= i) {
            int N = table.size();
            table.resize(i + 1);
            for (int j = N; j < int(table.size()); j++) {
                table[j] = mint(j).inv();
            }
        }
        return table[i];
    }
};

// ntt_friendlyå°‚ç”¨ 
template<typename mint>
struct fps : std::vector<mint> {
    using std::vector<mint>::vector;
    using self_t = fps<mint>;
    using mint_t = mint;

    fps(const std::vector<mint> &v) : std::vector<mint>::vector(v) {}

    int size() const {
        return std::vector<mint>::size();
    }

    bool empty() const {
        return size() == 0;
    }

    self_t inv(int K = -1) const {
        int N = size();
        assert(N != 0);
        if (K == -1) K = N;
        self_t g{(*this)[0].inv()}, tmp, tmp2;
        int step = 0;
        for (int k = 1; k < K; k *= 2, step++) {
            int sz = 4 * k;
            self_t tmp = g;
            self_t tmp2(this->begin(), this->begin() + std::min(N, 2 * k));
            tmp.resize(sz, 0);
            tmp2.resize(sz, 0);
            butterfly(tmp);
            butterfly(tmp2);
            for (int j = 0; j < sz; j++) tmp[j] = tmp[j] * tmp[j] * tmp2[j];
            butterfly_inv(tmp);
            sz = std::min(K, 2 * k);
            g.resize(sz, 0);
            for (int j = k; j < sz; j++) {
                g[j] = -tmp[j] * fps_base<mint>::ipow2(step + 2);
            }
        }
        return g;
    }

    self_t prefix(int K = -1) const {
        if (K == -1) K = size();
        K = std::min(K, size());
        return self_t(std::vector<mint>(this->begin(), this->begin() + K));
    }

    self_t suffix(int K = -1) const {
        int s = size();
        if (K == -1) K = s;
        K = std::min(K, s);
        return self_t(std::vector<mint>(this->begin() + s - K, this->end()));
    }

    self_t reverse(int K = -1) const {
        if (K == -1) K = size();
        K = std::min(K, size());
        self_t res(std::vector<mint>(this->begin(), this->begin() + K));
        std::reverse(res.begin(), res.end());
        return res;
    }

    self_t diff() const {
        int N = size();
        self_t res(std::max(0, N - 1));
        for (int i = 1; i < N; i++) {
            res[i - 1] = (*this)[i] * i;
        }
        return res;
    }

    self_t integral() const {
        int N = size();
        self_t res(N + 1);
        res[0] = 0;
        for (int i = 0; i < N; i++) {
            res[i + 1] = (*this)[i] * fps_base<mint>::inv_low(i + 1); 
        }
        return res;
    }

    self_t log(int K = -1) const {
        assert(!empty() && (*this)[0] == 1);
        if (K == -1) K = size();
        return (diff() * inv(K)).prefix(K - 1).integral();
    }

    self_t exp(int K = -1) const {
        assert(!empty() && (*this)[0] == 0);
        if (K == -1) K = size();
        self_t f{1}, g{1};
        int step = 0;
        self_t d = diff();
        for (int k = 1; k < K; k *= 2, step++) {
            mint i2 = fps_base<mint>::ipow2(step + 1);
            self_t fd = f.diff();
            self_t fpre = f;
            f.resize(2 * k, 0), g.resize(2 * k, 0);
            self_t tmp = g;
            butterfly(tmp);
            butterfly(f);
            for (int j = 0; j < 2 * k; j++) tmp[j] = tmp[j] * tmp[j] * f[j];
            butterfly_inv(tmp);
            for (int j = 0; j < k; j++) {
                g[j] = 2 * g[j] - tmp[j] * i2;
            }
            self_t q(d.begin(), d.begin() + std::min(d.size(), k - 1));
            self_t w = q;
            q.resize(2 * k, 0);
            butterfly(q);
            fd.resize(2 * k, 0);
            butterfly(fd);
            for (int j = 0; j < 2 * k; j++) q[j] = fd[j] - q[j] * f[j];
            w.resize(2 * k - 1, 0);
            self_t G = g;
            butterfly(G);
            for (int j = 0; j < 2 * k; j++) G[j] *= q[j];
            butterfly_inv(G);
            for (int j = k - 1; j < 2 * k - 1; j++) w[j] = G[j] * i2;
            w = w.integral();
            for (int j = 0; j < 2 * k; j++) {
                mint x = (j < size() ? (*this)[j] : 0);
                w[j] = x - w[j];
            }
            w[0]++;
            butterfly(w);
            for (int j = 0; j < 2 * k; j++) w[j] *= f[j];
            butterfly_inv(w);
            for (int j = 0; j < k; j++) f[j] = fpre[j];
            for (int j = k; j < 2 * k; j++) f[j] = w[j] * i2;
        }
        return f.prefix(K);
    }

    // f^b(x)
    self_t pow(long long b, int K = -1) const {
        assert(b >= 0);
        int N = size();
        if (K == -1) K = N;
        if (!b) return {1};
        int i = 0;
        while (i < N && (*this)[i] == 0) i++;
        long long min_deg = i;
        if (__builtin_mul_overflow(min_deg, b, &min_deg) || min_deg >= K || i >= N) {
            return self_t(K, 0);
        }
        mint iinv = (*this)[i].inv();
        self_t res(K - min_deg, 0);
        for (int j = i; j < std::min(N, (int)res.size() + i); j++) {
            res[j - i] = (*this)[j] * iinv;
        }
        res = res.log(-1);
        for (auto &x : res) x *= b;        
        res = res.exp(-1);
        mint ipow = (*this)[i].pow(b);
        for (auto &x : res) x *= ipow;
        self_t tmp(min_deg, 0);
        res.insert(res.begin(), tmp.begin(), tmp.end());
        return res;
    }

    std::pair<bool, self_t> sqrt(int K = -1) const {
        int N = size();
        if (K == -1) K = N;
        int d = 0;
        while (d < N && (*this)[d] == 0) d++;
        if (d == N || d / 2 >= K) return {true, self_t(K, 0)};
        if (d & 1) return {false, {}};
        auto [f, b0] = sqrt_mod<mint>((*this)[d]);
        if (!f) return {false, {}};
        self_t res{b0};
        mint i2 = fps_base<mint>::ipow2(1);
        K -= d / 2;
        for(int i = 1; i < K; i *= 2) {
            int M = std::min(K, 2 * i);
            res.resize(M, 0);
            self_t resi = res.inv(M) * self_t(this->begin() + d, this->begin() + std::min(N, d + M));
            for (int j = 0; j < M; j++) {
                res[j] = (res[j] + resi[j]) * i2;
            }
        }
        res.resize(K);
        self_t tmp(d / 2, 0);
        res.insert(res.begin(), tmp.begin(), tmp.end());
        return {true, res};
    }

    // q(x)g(x)+r(x) == thisã‚’æº€ãŸã™{q(x), r(x)} å¤šé …å¼ã®é™¤ç®—
    std::pair<self_t, self_t> division_polynomial(const self_t &g) const {
        int N = size(), M = g.size(), K = N - M + 1;
        assert(M != 0);
        if (N < M) return {{}, *this};
        self_t A(K), B(std::min(M, K), 0);
        for (int i = 0; i < K; i++) A[i] = (*this)[N - 1 - i];
        for (int i = 0; i < std::min(K, M); i++) B[i] = g[M - 1 - i];
        if (A.size() > K) A.resize(K);
        A *= B.inv(K);
        B.resize(K);
        for (int i = 0; i < K; i++) B[K - 1 - i] = A[i];
        A = B * g;
        A.resize(M - 1);
        for (int i = 0; i < M - 1; i++) A[i] = (*this)[i] - A[i];
        while (!A.empty() && A.back() == 0) A.pop_back();
        while (!B.empty() && B.back() == 0) B.pop_back();
        return {B, A};
    }

    // operator
    self_t operator - () const { self_t res(*this); for (int i = 0; i < size(); i++) res[i] = -res[i]; return res; }
    self_t operator += (const self_t &B) {
        if (size() < B.size()) this->resize(B.size(), 0);
        for (int i = 0; i < int(B.size()); i++) (*this)[i] += B[i];
        return *this;
    }
    self_t operator + (const self_t &B) const { self_t res(*this); return res += B; }
    self_t operator -= (const self_t &B) {
        if (size() < B.size()) this->resize(B.size(), 0);
        for (int i = 0; i < int(B.size()); i++) (*this)[i] -= B[i];
        return *this;
    }
    self_t operator - (const self_t &B) const { self_t res(*this); return res -= B; }
    self_t operator * (const self_t &B) const { return self_t(convolution_fft<mint>(*this, B)); }
    self_t operator *= (const self_t &B) { return *this = *this * B; }
    self_t operator /= (const self_t &B) { return *this *= B.inv(); }
    self_t operator / (const self_t &B) const { self_t res(*this); return res /= B; }
    self_t operator %= (const self_t &B) { return *this = division_polynomial(B).second; }
    self_t operator % (const self_t &B) const { self_t res(*this); return res %= B; }

    self_t operator += (const mint &B) {
        if (size() == 0) this->push_back(B);
        else (*this)[0] += B;
        return *this;
    }
    self_t operator + (const mint &B) const { self_t res(*this); return res += B; }
    self_t operator -= (const mint &B) {
        if (size() == 0) this->push_back(-B);
        else (*this)[0] -= B;
        return *this;
    }
    self_t operator - (const mint &B) const { self_t res(*this); return res -= B; }
    self_t operator * (const mint &B) const { self_t res(*this); return res *= B; }
    self_t operator *= (const mint &B) { for (int i = 0; i < size(); i++) (*this)[i] *= B; return *this; }
    self_t operator /= (const mint &B) { mint Bi = B.inv(); for (int i = 0; i < size(); i++) (*this)[i] *= Bi; return *this;  }
    self_t operator / (const mint &B) const { self_t res(*this); return res /= B; }
};

template<typename mint>
struct combination_mod {
  private:
    static int N;
    static std::vector<mint> F, FI, I;

  public:
    static bool built() { return !F.empty(); }

    static void clear() { N = 0; }

    // [0, N]ã‚’æ‰±ãˆã‚‹ã‚ˆã†ã«ã™ã‚‹
    // dynamic modint ç­‰ã§modã‚’å¤‰ãˆã¦å†ã³buildã™ã‚‹ã¨ãã¯clearã‚’å‘¼ã‚“ã§ãŠã
    // O(logMOD + å¢—ãˆãŸåˆ†)
    static void build(int _N) {
        _N++;
        assert(0 < _N && _N <= mint::mod());
        if (N >= _N) return;
        
        int preN = N;
        N = _N;
        F.resize(N);
        FI.resize(N);
        I.resize(N);
    
        F[0] = 1;
        for (int i = std::max(1, preN); i < N; i++) {
            F[i] = F[i - 1] * i;
        }
        FI[N - 1] = mint(F[N - 1]).inv();
        
        for (int i = N - 1; i >= std::max(1, preN); i--) {
            FI[i - 1] = FI[i] * i;
            I[i] = FI[i] * F[i - 1];
        }
    }

    static mint inv(int k) {
        return I[k];
    }

    using TypeMod = typename std::invoke_result<decltype(&mint::mod)>::type; // modintã®å†…éƒ¨çš„ãªæ•´æ•°åž‹
    
    static mint inv_large(TypeMod k) {
        if constexpr (std::is_same<TypeMod, int>::value) {
            long long res = 1;
            while (k >= N) {
                int q = -(mint::mod() / k);
                res *= q;
                res %= mint::mod();
                k = mint::mod() + q * k;
            }
            return mint(res) * I[k];
        } else {
            mint res = 1;
            while (k >= N) {
                TypeMod q = -(mint::mod() / k);
                res *= q;
                k = mint::mod() + q * k;
            }
            return res * I[k];
        }
    }

    static mint fac(int k) {
        return F[k];
    }

    static mint ifac(int k) {
        return FI[k];
    }

    static mint comb(int a, int b) {
        if (a < b || b < 0) return 0;
        return F[a] * FI[a - b] * FI[b];
    }

    static mint icomb(int a, int b) {
        assert(a >= b && b >= 0);
        return FI[a] * F[a - b] * F[b];
    }
    
    // O(b)
    static mint comb_small(int a, int b) {
        assert(b < mint::mod());
        if (a < b) return 0;
        mint res = 1;
        for (int i = 0; i < b; i++) res *= a - i;
        return res * FI[b];
    }

    // O(|b|) sum(b) = a
    static mint comb_multi(int a, const std::vector<int> &b) {
        mint res = 1;
        for (int r : b) {
            res *= comb(a, r);
            a -= r;
        }
        if (a == 0) return res;
        return 0;
    }

    static mint perm(int a, int b) {
        if (a < b || b < 0) return 0;
        return F[a] * FI[a - b];
    }

    static mint iperm(int a, int b) {
        assert(a >= b && b >= 0);
        return FI[a] * F[a - b];
    }

    // O(b)
    static mint perm_small(int a, int b) {
        assert(b < mint::mod());
        if (a < b) return 0;
        mint res = 1;
        for (int i = 0; i < b; i++) res *= a - i;
        return res;
    }
};

template<typename mint>
int combination_mod<mint>::N = 0;
template<typename mint>
std::vector<mint> combination_mod<mint>::F;
template<typename mint>
std::vector<mint> combination_mod<mint>::FI;
template<typename mint>
std::vector<mint> combination_mod<mint>::I;

// ä¸‹é™å†ª
// a0 + a1x + a2x(x-1) + .... anx(x-1)...(x-n+1)
template<typename mint>
struct falling_factorial : std::vector<mint> {
    using std::vector<mint>::vector;
    using self_t = falling_factorial<mint>;

    falling_factorial(const std::vector<mint> &v) : std::vector<mint>::vector(v) {}

    int size() const {
        return this->std::vector<mint>::size();
    }

    // f(0), f(1), .... f(N) ã‹ã‚‰Næ¬¡ã®falling factorialã‚’å¾©å…ƒ
    static self_t interpolation(const std::vector<mint> &y) {
        assert(!y.empty());
        int N = int(y.size()) - 1;
        combination_mod<mint>::build(N);
        std::vector<mint> X(N + 1), Y(N + 1);
        for (int i = 0; i <= N; i++) {
            X[i] = y[i] * combination_mod<mint>::ifac(i);
            Y[i] = (i % 2 == 0 ? 1 : -1) * combination_mod<mint>::ifac(i);
        }
        X = convolution_fft<mint>(X, Y);
        X.resize(N + 1);
        return self_t(X);
    }

    // f(0), f(1), .... f(K)ã‚’å¾—ã‚‹
    std::vector<mint> multipoint_evaluation(int K) const {
        combination_mod<mint>::build(K);
        std::vector<mint> y(K + 1);
        for (int i = 0; i <= K; i++) y[i] = combination_mod<mint>::ifac(i);
        auto res = convolution_fft<mint>(*this, y);
        res.resize(K + 1);
        for (int i = 0; i <= K; i++) res[i] *= combination_mod<mint>::fac(i);
        return res;
    }

    // f(x + c)
    self_t shift(mint c) const {
        int N = size() - 1;
        combination_mod<mint>::build(N);
        std::vector<mint> X(N + 1), Y(N + 1);
        mint C = 1;
        for (int i = 0; i <= N; i++) {
            X[i] = (*this)[i] * combination_mod<mint>::fac(i);
            Y[N - i] = C * combination_mod<mint>::ifac(i);
            C *= c - i;
        }
        X = convolution_fft<mint>(X, Y);
        self_t res(N + 1);
        for (int i = 0; i <= N; i++) {
            res[i] = X[N + i] * combination_mod<mint>::ifac(i);
        }
        return res;
    }
};

template<typename mint>
struct factorial {
    static std::vector<mint> bfac;
    static constexpr int K = 9;
    
    static void build() {
        bfac = {1, 3};
        for (int k = 1; k < K; k++) {
            auto f = falling_factorial<mint>::interpolation(bfac);
            f = f.shift(1 << k);
            auto tmp = f.multipoint_evaluation(3 << k);
            bfac.insert(bfac.end(), tmp.begin(), tmp.end());
            for (int j = 0; j < (2 << k); j++) {
                bfac[j] = bfac[2 * j] * bfac[2 * j + 1] * mint(2 * j + 1) * mint(1 << k);
            }
            bfac.resize(2 << k);
        }
        auto f = falling_factorial<mint>::interpolation(bfac);
        bfac = f.multipoint_evaluation((mint::mod() >> K) + 1);
        for (int i = 0; i < int(bfac.size()); i++) {
            bfac[i] *= (i + 1) << K;
            if (i) bfac[i] *= bfac[i - 1];
        }
    }
    
    // k!
    static mint fac(int k) {
        int block = k >> K;
        mint res = (block ? bfac[block - 1] : 1);
        for (int i = (block << K) + 1; i <= k; i++) res *= i;
        return res;
    }

    static mint perm(int n, int k) {
        if (k < 0 || n < k) return 0;
        return fac(n) / fac(k);
    }
    
    static mint comb(int n, int k) {
        if (k < 0 || n < k) return 0;
        return fac(n) / (fac(k) * fac(n - k));
    }
};
template<typename mint>
std::vector<mint> factorial<mint>::bfac;

int main() {
    factorial<modint998244353>::build(); 
    io_init();
    
    std::ofstream out("out.txt");
    
    if(out.is_open()) {
         for (int i = 1; i <= 10000000; i++) {
           out << factorial<modint998244353>::fac(i) << "," <<'\n';
         }
         out.close();
      } else {
           std::cerr << "Failed writing\n";
           return (1);
        }
    return EXIT_SUCCESS;
}
