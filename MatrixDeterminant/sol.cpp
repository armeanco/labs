#include <vector>
#include <cmath>
#include <utility>

long long determinant(std::vector< std::vector<long long> > m) {
  std::vector<std::vector<long double>> matrix(m.size(), std::vector<long double>(m.size()));
    for(std::size_t i = 0; i < m.size(); ++i) {
      for(std::size_t j = 0; j < m.size(); ++j) {
        matrix[i][j] = static_cast<long double>(m[i][j]);
      }
    }
    long double ans = 1.0;
    for(std::size_t i = 0; i < matrix.size(); ++i) {
      long double dif = 0.0;
      if(std::abs(matrix[i][i]) < 1e-12) {
        for(std::size_t s = i + 1; s < matrix.size(); ++s) {
          if(std::abs(matrix[s][i]) > 1e-12) {
            std::swap(matrix[i], matrix[s]);
            ans *= -1;
          }
        }
      }
      if(std::abs(matrix[i][i]) < 1e-12) return (0);
      for(std::size_t j = i + 1; j < matrix.size(); ++j) {
        dif = (matrix[j][i]) / (matrix[i][i]);
        for(std::size_t k = i; k < matrix[i].size(); ++k) {
          matrix[j][k] = matrix[j][k] - dif*matrix[i][k];
        }
      }
      ans *= matrix[i][i];
    }
    return std::llround(ans);
}
