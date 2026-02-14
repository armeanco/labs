def solve(lst):
    n = len(lst)
    
    hash = {}
    vis = {}
    
    sum = 0
    tmp = 0
    
    for i in range(n):
        if lst[i] not in hash:
            hash[lst[i]] = []
            vis[lst[i]] = 0
            tmp += (lst[i] * (n - i))
        
        hash[lst[i]].append(n - i)
        sum += lst[i]

    tmp -= lst[0]
    vis[lst[0]] += 1
    hash[lst[0]][0] -= 1
    
    ans = tmp
    for i in range(1, n - 1):
        vis[lst[i]] += 1
        if lst[i] != lst[i - 1]:
            tmp -= ((lst[i - 1] * hash[lst[i - 1]][vis[lst[i - 1]] - 1]) + lst[i])
            
            if hash[lst[i]][vis[lst[i]] - 1] > 0:
                hash[lst[i]][vis[lst[i]] - 1] -= 1
                hash[lst[i - 1]][vis[lst[i - 1]] - 1] -= 1
            
            if len(hash[lst[i - 1]]) > vis[lst[i - 1]]:
                tmp += (lst[i - 1] * hash[lst[i - 1]][vis[lst[i - 1]]])
                
            ans += tmp
            
        else:
            if hash[lst[i]][vis[lst[i]] - 1] > 0:
                hash[lst[i]][vis[lst[i]] - 1] -= 1
                
            tmp -= ((lst[i - 1] * hash[lst[i - 1]][vis[lst[i - 1]] - 1]) + lst[i])
            tmp += (lst[i] * (n - (i + 1)))
            ans += tmp

    return ans + sum
print(solve([]))
