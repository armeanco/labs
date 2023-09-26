#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <unordered_map>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <utility>
#include <numeric>
#include <map>
#include <climits>
#include <array>
#include <queue>
#include <cctype>

struct G {
   int x, y, keys;
   inline G(int a, int b, int c) noexcept : x(a), y(b), keys(c) { }

   bool operator <(const G& other) const {
      if( x != other.x ) return x < other.x;
      if( y != other.y ) return y < other.y;
      return keys < other.keys;
   }
};

int R, C;
std::vector<std::string> mp(101);
std::vector<std::vector<std::vector<int>>> visited(101, std::vector<std::vector<int>>(101, std::vector<int>(16)));
const int dx[4] = {0, 0, 1, -1}, dy[4] = {1, -1, 0, 0};

int sol(std::pair<int, int> &pos) {
	      int vmask = 0;
        std::queue<G> q;
        q.push({pos.first, pos.second, 0});
        visited[q.front().x][q.front().y][0] = 0;
        std::pair<int, int> v;
        while( !q.empty() ) {
           int u = q.front().x, t = q.front().y;
           int mask = q.front().keys;
           q.pop();
           for( int i = 0; i < 4; ++i ) {
              v = {u + dx[i], t + dy[i]};
              if( v.first >= 0 && v.first < R && v.second >= 0 && v.second < C && mp[v.first][v.second] != '#' ) {
                 if( mp[v.first][v.second] == 'X' ) return visited[u][t][mask] + 1;
                 if( visited[v.first][v.second][mask] == -1 ) {
                    vmask = mask;
                    if( std::islower(mp[v.first][v.second]) ) {
                       switch(mp[v.first][v.second]) {
                          case 'b': vmask |= 1; break;
                          case 'y': vmask |= 2; break;
                          case 'r': vmask |= 4; break;
                          case 'g': vmask |= 8; break;
                       }
                       q.push({v.first, v.second, vmask});
                       visited[v.first][v.second][mask] = visited[u][t][mask] + 1;
                       visited[v.first][v.second][vmask] = visited[v.first][v.second][mask];
                    }
                    else if( std::isupper(mp[v.first][v.second]) ) {
                       switch(mp[v.first][v.second]) {
                          case 'B': if( mask & 1 ) q.push({v.first, v.second, mask}); break;
                          case 'Y': if( mask & 2 ) q.push({v.first, v.second, mask}); break;
                          case 'R': if( mask & 4 ) q.push({v.first, v.second, mask}); break;
                          case 'G': if( mask & 8 ) q.push({v.first, v.second, mask}); break;
                       }
                       visited[v.first][v.second][mask] = visited[u][t][mask] + 1;
                    }
                    else {
                       q.push({v.first, v.second, mask});
                       visited[v.first][v.second][mask] = visited[u][t][mask] + 1;
                    }
                 }
              }
           }
        }
        return -1;
}

int main() {
  std::ios::sync_with_stdio(false);
	std::cin.tie(nullptr);
  std::cout.tie(nullptr);
  std::pair<int, int> pos;
  while( (std::cin >> R >> C) && ( R || C ) ) {
    for( int i = 0; i < R; ++i ) {
    std::cin >> mp[i];
    for( int j = 0; j < C; ++j ) { 
       if( mp[i][j] == '*' ) pos = {i, j};
          for( int k = 0; k < 16; ++k ) visited[i][j][k] = -1;
       }
    }
    int res = sol(pos);
    if( res == -1 ) std::cout << "The poor student is trapped!\n";
    else std::cout << "Escape possible in " << res << " steps.\n";
  }
	return EXIT_SUCCESS;
}
