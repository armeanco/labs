#include <string>
#include <vector>
#include <algorithm>

bool validateBase(const std::string& num, unsigned int base) {
  int mx = 0, i = 0;
  for(; i < (int)num.size(); ++i) mx = std::max(mx, (int)num[i] - '0');
  return (mx <= 9 && mx < (int)base) || (mx > 9 && mx - 6 <= (int)base) ? true : false;
}
