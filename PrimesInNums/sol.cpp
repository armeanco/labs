#include <vector>
#include <utility>
#include <string>

class PrimeDecomp
{
public:
    template<typename T>
    static std::vector<std::pair<T, int>> factor(T num) {
      std::vector<std::pair<T, int>> factors;
      for( T i = 2; i * i <= num; ++i ) {
        if( num % i == 0 ) {
          factors.push_back({i, 0});
          while( num % i == 0 ) {
            num /= i;
            factors.back().second++;
          }
        }
      }
      if( num > 1 ) factors.push_back({num, 1});
      return factors; 
    }
    static std::string factors(int lst) {
      std::vector<std::pair<int, int>> fact = factor(lst);
      std::string ans = "";
      for(auto &[a, b] : fact) ans += "(" + std::to_string(a) + (b > 1 ? "**" + std::to_string(b) : "") + ")";
      return ans;
    }
};
