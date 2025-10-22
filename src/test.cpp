#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <queue>
#include <stack>
#include <utility>
#include <cstring>
#include <cmath>
#include <limits>
// #include<sys/resource.h>
#include<cstdio>
using namespace std;
#define ll long long int
#define pb push_back
#define rb pop_back
#define ti tuple<int, int, int>
#define pii pair<int, int>
#define pli pair<ll, int>
#define pll pair<ll, ll>
#define mp make_pair
#define mt make_tuple
#define F first
#define S second

int main()
{
    // rlimit rlim;
    // getrlimit(RLIMIT_STACK, &rlim);
    // rlim.rlim_cur = rlim.rlim_max;
    // setrlimit(RLIMIT_STACK, &rlim);

    // freopen("input.txt", "r", stdin);
    // freopen("output.txt", "w", stdout);

    cin.tie(0)->sync_with_stdio(0);

    vector<pii> v;
    int cnt = 0;
    v.pb({1, cnt++});

    cout << v[0].F << " " << v[0].S << " " << cnt << endl;
}