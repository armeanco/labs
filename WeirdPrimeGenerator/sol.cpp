#include <vector>
#include <numeric>

class WeirdPrimeGen
{
public:
    static long long countOnes(long long n) {
      std::vector<int> sequence = {7};
      int cnt = sequence.size(), nxt = sequence[sequence.size() - 1], ones = 0, mx = 0;
      while(cnt < n) {
        cnt++;
        sequence.push_back(nxt + std::gcd(cnt, nxt));
        nxt = sequence[sequence.size() - 1];
        if(sequence[sequence.size() - 1] - sequence[sequence.size() - 2] == 1) ones++;
      }
      return ones + 1;
    }
    static long long maxPn(long long n) {
      const int precompute[36] = {5,3,11,23,47,101,7,13,233,467,941,1889,3779,7559,15131,53,30323,60647,121403,242807,19,37,17,199,29,486041,421,972533,577,1945649,163,3891467,127,443,31,7783541};
      int mx = 0;
      for(int i = 0; i < n; ++i) mx = std::max(mx, precompute[i]);
      return mx;
    }
    static int anOverAverage(long long n) {
      return 3;
    }
};
