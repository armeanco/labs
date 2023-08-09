#include <iostream>
#include <vector>
#include <utility>
#include <string>
#include <cmath>
#include <algorithm>
#include <set>

int evenSubarrays(std::vector<int> &array, int step) {
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
	int step;
	std::cin >> step;
	std::vector<int> array = {1, 2, 2, 1, 3, 3, 5, 5, 6};
	std::vector<int> array = {1, 2, 2, 2, 1, 1, 2};
	std::vector<int> array = {1, 2, 2, 1, 3, 3, 5, 5, 6, 1, 2, 2, 2, 1, 1, 2, 2, 6, 6};
	std::cout << evenSubarrays(array, step); 
	return EXIT_SUCCESS;
}
