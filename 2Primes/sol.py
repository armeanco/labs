def f(x):
    if x <= 1:
        return 0
    if x <= 3:
        return 1
    if x % 2 == 0 or x % 3 == 0:
        return 0
    k = 5
    while k * k <= x:
        if x % k == 0 or x % (k + 2) == 0:
            return 0
        k += 6
    return 1

def permutational_primes(n_max, k_perms):
    mp = [[0 for x in range(0)] for y in range(100000)]
    cnt = [0] * 100000
    frq = []
    low = 1e6
    ans = 0
    for x in range(13, n_max + 1):
        if f(x):
            l = [int(d) for d in str(x)]
            l.sort(reverse=True)
            q = int("".join(map(str, l)))
            cnt[q] += 1
            mp[q].append(x)
    for x in range(len(cnt)):
        if cnt[x] == k_perms + 1:
            low = mp[x][0]
            frq.append(low)
            ans += 1
    frq.sort()
    if ans == 0:
        return [0, 0, 0]
    return [ans, frq[0], frq[len(frq) - 1]]
