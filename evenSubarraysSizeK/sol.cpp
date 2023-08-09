#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <cmath>
#include <algorithm>
#include <set>

int evenSubarrays(std::vector<int> &array, int step) {
        if( step % 2 != 0 ) return 0;
	int cnt = 0, r = step, res = 0, l = 0, nxt = 0;
	std::vector<int> hash(100000);
	for( ; l++ < step; ) {
		hash[array[l - 1]]++;
		if( hash[array[l - 1]] % 2 != 0 ) cnt++; 
		if( hash[array[l - 1]] % 2 == 0 ) cnt--;
	}
	nxt = cnt;
	while( r < static_cast<int>(array.size()) ) {
		hash[array[r]]++;
		hash[array[r - step]]--;
		if( hash[array[r]] % 2 != 0 && array[r] != array[r - step] ) { cnt++; cnt == 0 ? res++ : res; }
		if( hash[array[r]] % 2 == 0 && array[r] != array[r - step] ) { cnt--; cnt == 0 ? res++ : res; }
		if( array[r] == array[r - step] ) { cnt == 0 ? res++ : res; }
		if( hash[array[r - step]] % 2 != 0 && array[r] != array[r - step] ) { cnt++; cnt == 0 ? res++ : res; }
		if( hash[array[r - step]] % 2 == 0 && array[r] != array[r - step] ) { cnt--; cnt == 0 ? res++ : res; }
		r++;
	}
	return nxt == 0 ? res + 1 : res;
}
int main() {
	int n, step;
	std::cin >> n >> step;
	std::vector<int> array(n);
	for( int i = 0; i < n; ++i ) std::cin >> array[i];
	std::cout << evenSubarrays(array, step); 
	return EXIT_SUCCESS;
}
