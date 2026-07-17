#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int c(const void *a, const void *b) {
    int _a = *(const int *)a;
    int _b = *(const int *)b;
    return (_a > _b) - (_a < _b);
}

int is_middle_decimal(long num) {
    if (num <= 0) return (num == 0) ? 1 : 0; 

    int digits[21];
    int cnt = 0;

    while (num > 0) {
        digits[cnt++] = num % 10;
        num /= 10;
    }
  
    qsort(digits, cnt, sizeof(int), c);
  
    long high = 0, low = 0;
    for (int i = cnt - 1; i >= 0; i--) high = high * 10 + digits[i];

    int nxt = 0;
    while (nxt < cnt && digits[nxt] == 0) {
        nxt++;
    }

    if (nxt < cnt) {
        low = digits[nxt];
        
        for (int i = 0; i < nxt; i++) low = low * 10;
        
        for (int i = nxt + 1; i < cnt; i++) low = low * 10 + digits[i];
    }
    
    if (!((low + high) % 2)) {
        return 1;
    }
    
    return 0;
}
