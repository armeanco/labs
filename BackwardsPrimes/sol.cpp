class BackWardsPrime
{
public:
  static std::string backwardsPrime(long long start, long long end) {
    auto f = [=](int x) -> bool {
      if(x <= 1) return false;
      if(x <= 2) return true;
      if(x % 2 == 0 || x % 3 == 0) return false;
      for(int k = 5; k * k <= x; k += 6) if(x % k == 0 || x % (k + 2) == 0) return false;
      return true;
    };
    std::string ans = "";
    for(int i = start; i <= end; ++i) {
      std::string t = "";
      if(f(i)) {
        t = std::to_string(i), std::reverse(t.begin(), t.end());
        if(std::stoi(t) != i && f(std::stoi(t))) ans += std::to_string(i) + " ";
      }
    }
    if(ans.size() > 0) ans.pop_back();
    return ans;
  }
};
