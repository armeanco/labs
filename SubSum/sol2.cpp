#include <iostream>
#include <vector>
#include <unordered_map>

long long solve(const std::vector<int>& arr) {
    int n = arr.size();
    if (n == 0) return 0;

    std::unordered_map<int, std::vector<int>> hash;
    std::unordered_map<int, int> vis;

    long long sum = 0ll, tmp = 0ll, ans = 0ll;

    for (int i = 0; i < n; ++i) {
        if (hash.find(arr[i]) == hash.end()) {
            hash[arr[i]] = std::vector<int>();
            vis[arr[i]] = 0;
            tmp += (1LL * arr[i] * (n - i));
        }
        hash[arr[i]].push_back(n - i);
        total_sum += arr[i];
    }

    tmp -= arr[0];
    vis[arr[0]] += 1;
    hash[arr[0]][0] -= 1;

    ans = tmp;

    for (int i = 1; i < n - 1; ++i) {
        vis[arr[i]] += 1;
        int current_val = arr[i];
        int prev_val = arr[i - 1];

        if (current_val != prev_val) {
            tmp -= ((1LL * prev_val * hash[prev_val][vis[prev_val] - 1]) + current_val);

            if (hash[current_val][vis[current_val] - 1] > 0) {
                hash[current_val][vis[current_val] - 1] -= 1;
                hash[prev_val][vis[prev_val] - 1] -= 1;
            }

            if (hash[prev_val].size() > (size_t)vis[prev_val]) {
                tmp += (1LL * prev_val * hash[prev_val][vis[prev_val]]);
            }
            ans += tmp;
        } 
        else {
            if (hash[current_val][vis[current_val] - 1] > 0) {
                hash[current_val][vis[current_val] - 1] -= 1;
            }

            tmp -= ((1LL * prev_val * hash[prev_val][vis[prev_val] - 1]) + current_val);
            tmp += (1LL * current_val * (n - (i + 1)));
            ans += tmp;
        }
    }

    return ans + sum;
}

int main() {
    int n;
    std::cin >> n;
    std::vector<int> arr(n);
    for(int i = 0; i < n; ++i) std::cin >> arr[i];
    std::cout << solve(arr);
    return EXIT_SUCCESS;
}
