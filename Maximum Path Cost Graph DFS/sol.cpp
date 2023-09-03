#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <unordered_map>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <numeric>
#include <map>
#include <climits>
#include <array>

int e[11][11];
int used[11];
        
void dfs(int v, int s, int n, int& ans) {
     used[v] = 1;
     ans = std::max(ans, s); // CHANGE FOR MIN COST PATH
     for( int i = 1; i <= n; ++i ) if( !used[i] && e[v][i] ) dfs( i, s + e[v][i], n, ans );
     used[v] = 0;
}

void sol() {
	std::ios::sync_with_stdio(false);
	std::cin.tie(nullptr);
        std::cout.tie(nullptr);
        int n, m, ans;
        std::cin >> n >> m;
        for( int i = 0; i < m; ++i ) {
            int a, b, c;
            std::cin >> a >> b >> c;
            e[a][b]=e[b][a]=c;
        }
        for( int i = 1; i <= n; ++i ) dfs( i, 0, n, ans );
        std::cout << ans << '\n';

}

int main() {
	sol();
	return EXIT_SUCCESS;
}
