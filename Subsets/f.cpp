//Number of subsets of a set equals to \binom(n)(k)=\frac(n!)(k!(n-k)!)
//where b(f(n)) = nth term in pascal triangle.
//0	1
//1	1   1
//2	1   2   1
//3	1   3   3   1
//4	1   4   6   4   1
//5	1   5   10   10   5   1
int main() {
  int32_t n;
  std::cin >> n;
  auto subsets = [](int32_t n) -> int64_t {
    return 2^n-1;
  };
  subsets(n);
}
//Now we can calculate the number of elements in each subset using the following formula
//Sum|B|(B|A)=n*2^n-1 where B is a subset of A
int main() {
  int32_t n;
  std::cin >> n;
  auto cnt = [](int32_t n) -> int64_t {
    return n*2^n-1;
  };
}
