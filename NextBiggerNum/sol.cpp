#include <vector>
#include <string>
#include <algorithm>

long nextBigger(long n) {
    int from = 0, to = 0;
    std::vector<int> prep, setik;
    for(int i = 0; i < static_cast<int>(std::to_string(n).size()); ++i) prep.push_back(std::to_string(n)[i] - '0');
    for(int i = 0; i < static_cast<int>(prep.size()); ++i) {
      for(int j = static_cast<int>(prep.size()); j > i; --j) {
        if(prep[i] < prep[j] && j < static_cast<int>(prep.size()) && j > 0) {
          from = i, to = j;
          break;
        }
      }
    }
    if(to == 0 && from == 0) return -1;
    int tmp = prep[from];
    prep[from] = prep[to], prep[to] = tmp;
    std::string ans = "";
    std::sort(prep.begin() + from + 1, prep.end());
    for(const auto &c : prep) ans += std::to_string(c);
    return std::stol(ans);
}
