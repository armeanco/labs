#include <string>
#include <vector>

struct Mat {
    long long mat[2][2];
    Mat() { mat[0][0] = mat[0][1] = mat[1][0] = mat[1][1] = 0; }
};

inline void mul(const Mat& A, const Mat& B, Mat& res, long long m) {
    long long r00 = (A.mat[0][0] * B.mat[0][0] + A.mat[0][1] * B.mat[1][0]) % m;
    long long r01 = (A.mat[0][0] * B.mat[0][1] + A.mat[0][1] * B.mat[1][1]) % m;
    long long r10 = (A.mat[1][0] * B.mat[0][0] + A.mat[1][1] * B.mat[1][0]) % m;
    long long r11 = (A.mat[1][0] * B.mat[0][1] + A.mat[1][1] * B.mat[1][1]) % m;
    res.mat[0][0] = r00; res.mat[0][1] = r01;
    res.mat[1][0] = r10; res.mat[1][1] = r11;
}

Mat powersOfT[10];

void precompute(long long m) {
    Mat T;
    T.mat[0][0] = 1; T.mat[0][1] = 1;
    T.mat[1][0] = 1; T.mat[1][1] = 0;
    
    powersOfT[0].mat[0][0] = 1; powersOfT[0].mat[1][1] = 1;
    for(int i = 1; i < 10; i++) {
        mul(powersOfT[i-1], T, powersOfT[i], m);
    }
}

long long fib(std::string n, int p) {
    precompute(p);

    Mat ans;
    ans.mat[0][0] = 1; ans.mat[1][1] = 1;

    for (char c : n) {
        int dig = c - '0';
      
        Mat a2, a4, a8;
        mul(ans, ans, a2, p);
        mul(a2, a2, a4, p);
        mul(a4, a4, a8, p);
        mul(a8, a2, ans, p); 

        if (dig > 0) {
            Mat temp;
            mul(ans, powersOfT[dig], temp, p);
            ans = temp;
        }
    }
    return ans.mat[0][1];
}
