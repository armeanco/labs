#include <iostream>
#include <vector>

int vis[10000000];

long long sol(std::vector<int> lst) {

   std::vector<std::vector<int>> hash(10000000, std::vector<int>());
   int n = static_cast<int>(lst.size());
   
   long long sum = 0ll, tmp = 0ll, ans = 0ll;

   for(int i = 0; i < n; ++i) {
      hash[lst[i]].push_back(n - i);
      sum += lst[i];
      if(hash[lst[i]].size() < 2) tmp += (lst[i] * static_cast<long long>(n - i));
   }

   tmp -= lst[0];
   vis[lst[0]]++;
   hash[lst[0]][0] -= 1;
   
   for(int i = 1; i < n - 1; ++i) {
      vis[lst[i]]++;
      if(lst[i] != lst[i - 1]) {
        tmp -= ((lst[i - 1] * (hash[lst[i - 1]][vis[lst[i - 1]] - 1])) + lst[i]);
        if(hash[lst[i]][vis[lst[i]] - 1] > 0) { hash[lst[i]][vis[lst[i]] - 1]--; hash[lst[i - 1]][vis[lst[i - 1]] - 1]--; }
        if(hash[lst[i - 1]].size() > vis[lst[i - 1]] && i - 1 > 0) tmp += (lst[i - 1] * (hash[lst[i - 1]][vis[lst[i - 1]]]));
        ans += tmp;
      }
      else if(lst[i] == lst[i - 1]) {
        if(hash[lst[i]][vis[lst[i]] - 1] > 0) hash[lst[i]][vis[lst[i]] - 1]--;
        tmp -= ((lst[i - 1] * (hash[lst[i - 1]][vis[lst[i - 1]] - 1])) + lst[i]);
        tmp += (lst[i] * (n - (i + 1)));
        ans += tmp;
      }
   }

   return (ans + sum);
}

int main() {
  
  int n;
  std::cin >> n;

  std::vector<int> lst(n);

  for(int i = 0; i < n; ++i) {
    std::cin >> lst[i];
  }

  std::cout << sol(lst);

  return EXIT_SUCCESS;
}
