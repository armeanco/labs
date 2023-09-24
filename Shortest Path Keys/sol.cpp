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

struct G {
   int x, y, keys;
   G(int a, int b, int c) noexcept : x( a ), y( b ), keys( c ) { } 

   bool operator<(const G& other) const {
      if( x != other.x ) return x < other.x;
      if( y != other.y ) return y < other.y;
      return keys < other.keys;
   }
};

void sol(std::vector<std::string> &grid) {
  std::pair<int, int> start;
  int keyCnt = 0, m = static_cast<int>( grid.size() ), n = static_cast<int>( grid[0].size() );
  if( m == 0 || n == 0 ) return;
        
  for( int i = 0; i < m; ++i ) {
     for( int j = 0; j < n; ++j ) {
        if( grid[i][j] == '@' ) start = {i, j};
        else if( grid[i][j] >= 'a' && grid[i][j] <= 'z' ) keyCnt++;
     }
  }

  if( keyCnt == 0 ) return;
  std::queue<G> q;
  std::set<G> visited;
  q.push( {start.first, start.second, 0} );
  int steps = 0, allKeys = ( 1 << keyCnt ) - 1;
        
  while( !q.empty() ) {
     int size = q.size();
     for( int i = 0; i < size; ++i ) {
        G s = q.front(); q.pop();
        if( s.x < 0 || s.x >= m || s.y < 0 || s.y >= n || grid[s.x][s.y] == '#' ) continue;
        if (visited.count( s ) > 0 ) continue;
        visited.insert( s );
              
        if( grid[s.x][s.y] >= 'a' && grid[s.x][s.y] <= 'z' ) s.keys = ( ( s.keys ) | ( 1 << ( grid[s.x][s.y] - 'a' ) ) );
        else if( grid[s.x][s.y] >= 'A' && grid[s.x][s.y] <= 'Z' ) if( ( (s.keys) >> ( ( grid[s.x][s.y] - 'A' ) ) & 1) != 1 ) continue;
        visited.insert( s );

        if( s.keys == allKeys ) { std::cout << steps; return; }

        q.push( {s.x - 1, s.y, s.keys} );
        q.push( {s.x + 1, s.y, s.keys} );
        q.push( {s.x, s.y - 1, s.keys} );
        q.push( {s.x, s.y + 1, s.keys} );
     }
     steps++;
  }
  std::cout << "-1" << '\n';
}

int main() {
  std::ios::sync_with_stdio(false);
	std::cin.tie(nullptr);
  std::cout.tie(nullptr);
  int N;
  std::cin >> N;
  std::vector<std::string> grid( N );
  for( int i = 0; i < N; ++i ) std::cin >> grid[i];
	sol(grid);
	return EXIT_SUCCESS;
}
