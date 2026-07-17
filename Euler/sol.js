function totient(n) {
  if(n == 1) return 1;
  if(!Number.isInteger(n) || n <= 0) return 0;
  const f = (x) => {
    let res = n, p = 2;
    while(p * p <= x) {
        if(x % p == 0) {
            while(x % p == 0) x /= p;
            res -= res / p;
        }
        p += 1;
    }
    if(x > 1) res -= res / x;
    return res;
  };
  return f(n);
}
