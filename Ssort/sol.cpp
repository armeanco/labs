#include <vector>

void stalinSort(std::vector<int>& arr) {
  if(arr.size() > 0) {
    int mx = arr[0], i = 1;
    std::vector<int> cnt{arr[0]};
    for(; i < static_cast<int>(arr.size()); ++i) {
      if(mx <= arr[i]) { mx = arr[i]; cnt.push_back(mx); }
    }
    arr = cnt;
  }
}
