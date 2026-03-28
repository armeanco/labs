#include <vector>
#include <array>
#include <algorithm>

struct G {
public:
  unsigned int count = 0, low = 0;
};

std::array<unsigned int, 3> p_primes(const unsigned int n_max, const unsigned int k_perms) {
    auto f = [&](unsigned int x) -> bool {
        if (x <= 1) return false;
        if (x <= 3) return true;
        if (x % 2 == 0 || x % 3 == 0) return false;
        for (unsigned int i = 5; i * i <= x; i += 6) {
            if (x % i == 0 || x % (i + 2) == 0) return false;
        }
        return true;
    };
    
    auto precompute = [=](unsigned int n) -> unsigned int {
      int digits[10] = {0};
      while (n > 0) {
        digits[n % 10]++;
        n /= 10;
      }
      unsigned int d = 0;
      for (int i = 9; i >= 0; --i) {
        for (int j = 0; j < digits[i]; ++j) {
          d = d * 10 + i;
        }
      }
      return d;
    };
    
    std::vector<G> frequency(1000000);
    std::vector<unsigned int> left;
    
    std::array<unsigned int, 3> ans = {0, 1000000, 0};

    for (unsigned int i = 13; i <= n_max; ++i) {
        if (f(i)) {
            unsigned int d = precompute(i);
            if (frequency[d].count == 0) {
                frequency[d].low = i;
            }
            frequency[d].count++;
            if (frequency[d].count == k_perms + 1) {
                left.push_back(d);
            }
        }
    }
  
    for (unsigned int k : left) {
        if (frequency[k].count == k_perms + 1) {
            ans[0]++, ans[1] = std::min(ans[1], frequency[k].low), ans[2] = std::max(ans[2], frequency[k].low);
        }
    }

    return ans[0] == 0 ? std::array<unsigned int, 3>{0, 0, 0} : ans;
}
