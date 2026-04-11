//Optimal linear
#include <string>
#include <vector>
#include <cstdint>

std::uint64_t vowelRecognition(const std::string &str) {
    int cnt = 0;
    std::vector<int> frequency, precompute;
    long long carry = 0ll, range = 0ll;
    for(std::size_t i = 0; i < str.size(); ++i) {
      if(str[i] == 'A' || str[i] == 'E' || str[i] == 'U' || str[i] == 'I' || str[i] == 'O' ||
         str[i] == 'a' || str[i] == 'e' || str[i] == 'u' || str[i] == 'i' || str[i] == 'o') {
           cnt++;
         }
      if(cnt > 0) frequency.push_back(cnt);
    }
    for(int i = 0; i < static_cast<int>(frequency.size()); ++i) {
      int nxt = 1;
      if(frequency[i] == frequency[i + 1]) {
        int cur = frequency[i];
        while(frequency[i] == frequency[i + 1] && i + 1 < static_cast<int>(frequency.size())) {
          nxt++;
          i++;
          cur += frequency[i];
        }
        range += nxt;
        carry += cur;
        precompute.push_back(nxt);
      }
      else {
        precompute.push_back(1);
        carry += frequency[i];
        range += 1;
      }
    }
    long long res = 0ll, prev = 0ll;
    int already_taken = 0;
    res = carry;
    for(int i = 0; i < static_cast<int>(str.size()); ++i) {
      if(str[i] == 'A' || str[i] == 'E' || str[i] == 'U' || str[i] == 'I' || str[i] == 'O' ||
         str[i] == 'a' || str[i] == 'e' || str[i] == 'u' || str[i] == 'i' || str[i] == 'o') {
           res -= range;
           range -= precompute[already_taken];
           prev += res;
           already_taken++;
         }
         else {
           prev += res;
         }
    }
    return prev + carry;
}
