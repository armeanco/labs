#include <iostream>
#include <cmath>

#define MOD 999997

int main() {
    long long n, t, ans;
    std::cin >> n;
    auto g = [&]() -> long long {
      for(long long i = 1; i*i <= n; ++i) ans += (n / i) % MOD;
      t = n / (std::sqrt(n) + 1);
      for(long long i = 1; i <= t; ++i) ans += i * ((n / i) - (n / (i + 1))) % MOD;
      return ans;
    };
    g();
    std::cout << ans << '\n';
    return EXIT_SUCCESS;
}
