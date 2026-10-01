#include <iostream>
#include <string>
#include <vector>

int main() {
    long long a, b;
    std::cin >> a >> b;
    if(a == b) {
        std::cout << a << '\n';
        return 0;
    }
    int sum_a = 0, sum_b = 0, nxt = 0;
    int x = static_cast<int>(std::to_string(a).size()), y = static_cast<int>(std::to_string(b).size());
    std::vector<long long> hash(10000);
    for(int i = 0; i < (x < y ? y : x); ++i) {
        if(first > 0) {
            sum_a += first % 10;
            first /= 10;
        }
        if(second > 0) {
            sum_b += second % 10;
            second /= 10;
        }
    }
    for(int i = 0; i < y; ++i) {
        std::string left, right, tmp;
        tmp = (std::to_string(b)[i] - '0' - 1) >= 0 ? std::to_string(std::to_string(b)[i] - '0' - 1) : "x";
        long long rekt = 0;
        if(tmp != "x") {
            for(int j = 0; j < i; ++j) {
                left += std::to_string(b)[j];
            }
            for(int k = i + 1; k < y; ++k) {
                right += "9";
            }
        }
        if(tmp != "x") {
            left += tmp + right;
            long long qeq = std::stoll(left);
            if(qeq >= a) {
                long long f = qeq;
                while(qeq > 0) {
                    rekt += qeq % 10;
                    qeq /= 10;
                }
                hash[rekt] = f;
            }
            nxt = std::max(nxt, rekt);
        }
    }
    hash[sum_a] = a;
    hash[sum_b] = b;
    std::cout << (nxt > std::max(sum_a, sum_b) ? hash[nxt] : hash[std::max(sum_a, sum_b)]) << '\n';
    return 0;
}
