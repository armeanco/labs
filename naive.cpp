#include <iostream>
#include <array>
#include <cstdlib>
#include <random>

constexpr int MOD =  2099;
constexpr int SIZE = 450;
constexpr int N =    2e5+5;
constexpr int M =    1e5+5;
constexpr int B =    450;

int las[M];
int bel[N];
int L[B];
int R[B];
unsigned long long val[M];
unsigned long long seq[N];
unsigned long long tag[B];
std::mt19937_64 Rand(19511016);

class HashMap {
public:
	int N;
	int Time;               
	int Ticks               [MOD  + 2];
	int Hop                 [MOD  + 2];
	int next                [SIZE + 2];
	int count               [SIZE + 2];
	unsigned long long read [SIZE + 2];

	void A(int x, unsigned long long y, int z);

	void Free();

	void Set(int x);

	int Get(int x);

	int Queue(unsigned long long x);

	void Put(unsigned long long x);
} cnt[B];

void HashMap::A(int x, unsigned long long y, int z) {
	read[++N] = y;
	next[N] = Hop[x];
	count[N] = z;
	Hop[x] = N;
}

void HashMap::Free() {
	++Time;
	N = 0;
}

void HashMap::Set(int x) {
	HashMap::Free();
	HashMap::Get(0);
	HashMap::A(0, 0, x);
}

int HashMap::Get(int x) {
	return Ticks[x] == Time ? Hop[x] : (Ticks[x] = Time, Hop[x] = 0);
}

int HashMap::Queue(unsigned long long x) {
	for(int i = HashMap::Get(x % MOD); i; i = next[i]) {
			if(read[i] == x) {
				return count[i];
			}
	}

	return EXIT_SUCCESS;
}

void HashMap::Put(unsigned long long x) {
	for(int i = HashMap::Get(x % MOD); i; i = next[i]) {
		if(read[i] == x) {
			return void(++count[i]);
		}
	}
	HashMap::A((x % MOD), x, 1);
}

void Update(int x, unsigned long long y) {

	if(!x) {
		return;
	}

	int q = bel[x];
	for(int i = 1; i < q; ++i) {
		tag[i] ^= y;
	}

	cnt[q].Free();

	for(int i = L[q]; i <= x; ++i) {
		cnt[q].Put(seq[i]^=y);
	}

	for(int i = R[q]; i > x; --i) {
		cnt[q].Put(seq[i]);
	}

}

int Query(int x) {

	int q = bel[x];
	int ans = 0;
	unsigned long long v = tag[q];
	for(int i = 1; i < q; ++i) {
		ans += cnt[i].Queue(tag[i]);
	}
	
	for(int i = L[q]; i <= x; ++i) {
		ans += (seq[i] == v);
	}
	return ans;
}
	
template<std::size_t S>
long long naive(std::array<int, S> input) {

	int n = input.size();
	for(int i = 1; i <= n; ++i) {
		seq[i] = 0;
		las[input[i - 1]] = 0;
	}

	for(int i = 1; i <= n; ++i) {
		if(!val[input[i - 1]]) {
			val[input[i - 1]] = Rand();
		}
	}
	

	for(int i = 1; i <= n; ++i) {
		bel[i] = (i - 1)/SIZE + 1;
	}

	int tot = bel[n];

	for(int i = 1; i <= tot; ++i) {
		L[i] = (i - 1)*SIZE+1;
		R[i] = i *SIZE;
	}

	R[tot] = std::min(R[tot], n);

	for(int i = 1; i <= tot; ++i) {
		tag[i] = 0;
		cnt[i].Set(R[i] - L[i] + 1);
	}

	long long ans = 0;
	for(int i = 1; i <= n; ++i) {
		Update(las[input[i - 1]], val[input[i - 1]]);
		ans += Query(i);
		las[input[i - 1]] = i;
	}
	return ans;
}

int RUN__TEST__FIRST() {

	std::random_device dev;
	std::mt19937_64 rng(dev());
	std::uniform_int_distribution<std::mt19937_64::result_type> dist(1, 100);

	std::array<int, 200000> test__lib__drive;

	for(int i = 0; i < 200000; ++i) {
		test__lib__drive[i] = dist(rng);
	}

	return naive(test__lib__drive);
}

int RUN__TEST__SECOND() {

	return EXIT_SUCCESS;
}


int main(int argc, char** argv) {
	
	std::array<int, 4>  test__lib__one {2, 2, 2, 3};
	std::array<int, 15> test__lib__two {2, 5, 2, 3, 6, 7, 8, 23, 23, 13, 65, 31, 3, 4, 3};
	std::cout << naive(test__lib__one) << '\n';
	std::cout << naive(test__lib__two);

	std::cout << RUN__TEST__FIRST() << '\n';
	std::cout << RUN__TEST__FIRST() << '\n';

	return EXIT_SUCCESS;
}
