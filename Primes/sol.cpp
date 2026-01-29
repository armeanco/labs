//https://judge.yosupo.jp/problem/counting_primes
//Compute P(N) for N = 10^13 under ~3000ms
#include <stdio.h>
#include <math.h>
#include <algorithm>
#include <vector>
typedef long long ll;
inline ll piSieve(const ll N){
    register const ll n = N;
    if(n<=1) return 0LL;
    if(n==2) return 1LL;
    register const int lim=std::sqrt(n);
    int vsz=(lim+1)>>1;
    std::vector<int> smalls(vsz);
    for(int cx=0;cx<vsz;++cx) smalls[cx]=cx;
    std::vector<int> roughs(vsz);
    for(int cx=0;cx<vsz;++cx) roughs[cx]=(cx<<1|1);
    std::vector<ll> larges(vsz);
    for(int cx=0;cx<vsz;++cx) larges[cx]=(n/(cx<<1|1)-1)>>1;
    std::vector<bool> skips(lim+1);
    int pCnt=0;
    for(register int p=3;p<=lim;p+=2){
        if(skips[p]) continue;
        int p2=p*p;
        if(1LL*p2*p2>n) break;
        skips[p]=true;
        for(int cx=p2;cx<=lim;cx+=(p<<1))
            skips[cx]=true;
        int ns=0;
        for(int cx=0;cx<vsz;++cx){
            int cur=roughs[cx];
            if(skips[cur]) continue;
            ll d=1LL*cur*p;
            larges[ns]=larges[cx]-(d<=lim?larges[smalls[d>>1]-pCnt]
                                    :smalls[(ll((double(n)/d)-1))>>1])+pCnt;
            roughs[ns++]=cur;
        }
        vsz=ns;
        for(int cx=(lim-1)>>1,cy=((lim/p)-1)|1;cy>=p;cy-=2){
            int cur=smalls[cy>>1]-pCnt;
            for(int cz=(cy*p)>>1;cz<=cx;--cx)
                smalls[cx]-=cur;
        }
        ++pCnt;
    }
    larges[0]+=1LL*(vsz+((pCnt-1)<<1))*(vsz-1)>>1;
    for(int cx=1;cx<vsz;++cx) larges[0]-=larges[cx];
    for(int cx=1;cx<vsz;++cx){
        int q=roughs[cx];
        ll m=n/q;
        int e=smalls[((m/q)-1)>>1]-pCnt;
        if(e<cx+1) break;
        ll t=0;
        for(int cy=cx+1;cy<=e;++cy)
            t+=smalls[ll((double(m)/roughs[cy])-1)>>1];
        larges[0]+=t-1LL*(e-cx)*(pCnt+cx-1);
    }
    return larges[0]+1;
}
int main() {
    register ll n;
    ll k;
    scanf("%lld", &n);
    printf("%lld\n", piSieve(n));
    return 0;
}
