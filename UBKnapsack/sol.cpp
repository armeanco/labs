#include <iostream>
#include <vector>

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);

    int N, K;
    std::cin >> N >> K;
    std::vector<std::vector<bool>> table(N, std::vector<bool>(K + 1, 0));
    table[N - 1][0] = 1;
    for(int i = N; i <= K; i += N) {
        table[N - 1][i] = 1;
    }
    for(int i = N - 2; i >= 0; --i) {
        for(int j = 0; j <= K; ++j) {
            table[i][j] = table[i][j] || table[i + 1][j];
            if(j >= (i + 1)) {
                table[i][j] = table[i][j] || table[i][j - (i + 1)];
            }
        }
    }
    for(const auto &c : table) {
        for(const auto &p : c) {
            std::cout << p << " ";
        }
        std::cout << "\n";
    }
    return 0;
}
