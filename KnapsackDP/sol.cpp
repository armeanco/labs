#include <iostream>
#include <vector>
#include <algorithm>

int knap(int v, const std::vector<int> &choose, const std::vector<int> &sequence) {
    std::vector<std::vector<int>> dp(static_cast<int>(choose.size()) + 1, std::vector<int>(v + 1, 0));
    
    for (int i = 1; i <= static_cast<int>(choose.size()); ++i) {
        for (int w = 1; w <= v; ++w) {
            if (choose[i - 1] <= w)
                dp[i][w] = std::max(sequence[i - 1] + dp[i - 1][w - choose[i - 1]], dp[i - 1][w]);
            else
                dp[i][w] = dp[i - 1][w];
        }
    }
    return dp[static_cast<int>(choose.size())][v];
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    int n, v;
    std::cin >> n >> v;

    std::vector<int> sequence(n), choose(n);
    for(int i = 0; i < n; ++i) {
        std::cin >> sequence[i];
        choose[i] = i + 1;
    }

    int ans = knap(v, choose, sequence);
    std::cout << ans << '\n';
    return 0;
}
