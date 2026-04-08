#include <string>
#include <vector>
#include <utility>
#include <cstdint>

//Sub-optimal solution O(N*V)
std::uint64_t vowelRecognition(const std::string &str) {
    std::uint64_t cnt = 0, ans = 0, def = 0;
    std::vector<int> frequency, sequence, preprocess = {1};
    std::vector<std::pair<int, int>> compressed;
    for(std::size_t i = 0; i < str.size(); ++i) {
      if(str[i] == 'A' || str[i] == 'E' || str[i] == 'U' || str[i] == 'I' || str[i] == 'O' ||
         str[i] == 'a' || str[i] == 'e' || str[i] == 'u' || str[i] == 'i' || str[i] == 'o') {
           cnt++;
           sequence.push_back(i + 1);
         }
      if(cnt > 0) frequency.push_back(cnt);
    }
    for(int i = 0; i < frequency.size(); ++i) {
      std::uint64_t nxt = 1;
      if(frequency[i] == frequency[i + 1]) {
        int cur = frequency[i];
        while(frequency[i] == frequency[i + 1] && i + 1 < frequency.size()) {
          nxt++;
          i++;
          cur += frequency[i];
        }
        compressed.push_back(std::make_pair(nxt, cur));
        ans += cur;
      }
      else {
        compressed.push_back(std::make_pair(1, frequency[i]));
        ans += frequency[i];
      }
    }
    def = ans;
    for(int i = 1; i < sequence.size(); ++i) {
      preprocess.push_back(sequence[i] - sequence[i - 1]);
    }
    for(int i = 0; i < preprocess.size(); ++i) {
      std::uint64_t sum = 0;
      for(int j = i; j < compressed.size(); ++j) {
        compressed[j].second -= compressed[j].first;
        sum += compressed[j].second >= 0 ? compressed[j].second : 0;
      }
      sum *= preprocess[i + 1] * 1ll;
      ans += sum;
      sum = 0;
    }
    if(sequence.size() && sequence[0] > 1) ans += def * (sequence[0] - 1ll);
    return ans;
}
