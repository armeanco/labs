#include <vector>

unsigned int findBall(Scales s) {
  auto q = [&](std::initializer_list<int> a, std::initializer_list<int> b) { 
      return s.getWeight(a, b); 
  };
  int i = q({0, 3, 6}, {1, 4, 7}), j = q({0, 1, 2}, {3, 4, 5});
  return (((2 - i) % 3 + 3) % 3) + 3 * (((2 - j) % 3 + 3) % 3);
}
