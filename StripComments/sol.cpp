#include <string>
#include <unordered_set>
#include <regex>

std::string stripComments(const std::string& str, const std::unordered_set<char>& markers) {
  int cnt = 0;
  std::string s = str;
  for(const auto &c : markers) {
    for(int i = 0; i < static_cast<int>(s.size()); ++i) {
      if(s[i] == c) {
        s[i] = '9';
        if(i - 1 > 0) s[i - 1] = '9';
        for(int j = i + 1; s[j] != '\n' && j < static_cast<int>(s.size()); ++j, ++i) s[j] = '9';
      }
      else if(s[i] == '\n' && s[i - 1] == ' ' && i - 1 > 0) s[i - 1] = '9';
      else if(s[i] == ' ') cnt++;
    }
    if(cnt == static_cast<int>(str.size())) return "";
  }
  s = std::regex_replace(s, std::regex("\\d+"), "");
  if(s[s.size() - 1] == ' ' && s.size() > 1) s.pop_back();
  return s;
}
