#include <math.h>
#include <stdio.h>

unsigned long long sum_1_to_n (unsigned n)
{
    return n * (n + 1ull) / 2ull;
}

void cantor_pairing_function(char *output, unsigned n)
{
// solve 1 + 2 + 3 + ... + diagonal > n (quadratic equation)
    unsigned diagonal = ceil( (-1 + sqrt(8 * n + 1) ) / 2);
    unsigned position = n - sum_1_to_n(diagonal - 1);
// odd diagonals: d/1 to 1/d; even diagonals: 1/d to d/1
    unsigned numerator = (diagonal % 2) ? (diagonal - position + 1) : position;
// on diagonal d, (numerator + denominator) == d + 1
    unsigned denominator = diagonal + 1 - numerator;
    sprintf(output, "%u/%u", numerator, denominator);
}
