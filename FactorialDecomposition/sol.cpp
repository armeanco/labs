#include <string>
#include <vector>

std::string decomp(int n) {
  std::string ans = "";
  auto f = [=](int n) -> std::vector<int> {
    if (n < 2) return {};
    int size = (n - 1) / 2;
    std::vector<bool> ok(size + 1, true);
    std::vector<int> sequence;
    for (int i = 1; i * i <= n; ++i) {
        if (ok[i]) {
            int p = 2 * i + 1;
            for (int j = (p * p - 1) / 2; j <= size; j += p) ok[j] = false;
        }
    }
    sequence.push_back(2);
    for (int i = 1; i <= size; ++i) if (ok[i]) sequence.push_back(2 * i + 1);
    return sequence;
  };
  std::vector<int> precompute = f(n);
  int cnt = 0, nxt = 0;
  for(int i = 0; i < static_cast<int>(precompute.size()) - 1; ++i) {
    for(int j = n; j >= 1; j /= precompute[i], cnt++, nxt += j) {}
    if(nxt == 1) ans += std::to_string(precompute[i]) + " * ";
    else ans += std::to_string(precompute[i]) + "^" + std::to_string(nxt) + " * ";
    cnt = 0, nxt = 0;
  }
  ans += std::to_string(precompute[precompute.size() - 1]);
  return ans;
}
