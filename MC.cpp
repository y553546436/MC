#include <cstdio>
#include <algorithm>
#include <cstring>
#include <vector>
#include <utility>
#include <ctime>
#include <list>
using namespace std;
const int N = 755;
struct edge {
    int lastuncover,u,v;
}; 
list<edge> L, UL;

bool G[N][N], C[N], _C[N], covered[N][N];
int Csize, Lsize;
int n, m, dscore[N], w[N][N], lastuncover[N][N];

int step;

void init() {
    L.clear();
    UL.clear();
    memset(G,0,sizeof G);
    memset(C,0,sizeof C);
    memset(_C, 0, sizeof _C);
    memset(covered, 0, sizeof covered);
    Csize=Lsize=0;
    memset(lastuncover,0,sizeof lastuncover);
}

void addtoC(int u) {
    C[u] = 1; Csize++;
    dscore[u] = -dscore[u];
    for (int i = 0; i < n; ++i)
        if (G[u][i]) {
            if (C[i]) dscore[i] += w[u][i];
            else {
                dscore[i] -= w[u][i];
                Lsize --;
                covered[u][i] = true;
            }
        }
}
void removefromC(int u) {
    C[u] = 0; Csize--;
    dscore[u] = -dscore[u];
    for (int i = 0; i < n; ++i)
        if(G[u][i]) {
            if (C[i]) dscore[i] -= w[u][i];
            else {
                dscore[i] += w[u][i];
                Lsize ++;
                covered[u][i] = false;
                lastuncover[u][i] = step;
                L.push_back((edge){step,u,i});
                UL.push_back((edge){step,u,i});
            }
        }
}

void remove_some(int ub, int delta) {
    while (Csize > 1 && Csize > ub - delta) {
        //printf("Csize=%d ub=%d\n", Csize, ub);
        int index = rand() % Csize;
        int cnt = 0;
        for (int i = 0; i < n; ++i)
            if (C[i] && (cnt++) == index) {
                removefromC(i);
                break;
            }
    }
}

bool ChooseExchangePair(int &last_add, int &last_remove) {
    int p = -1, q = -1;
    auto it = L.begin();
    while(it != L.end()) {
        int u = it->u, v = it->v;
        if (covered[u][v] || it->lastuncover != lastuncover[u][v])
            it = L.erase(it);
        else {
            for (int i = 0; i < n; ++i)
                if (C[i]) {
                    int score_u = dscore[i] + dscore[u] + w[i][u];
                    int score_v = dscore[i] + dscore[v] + w[i][v];
                    if (score_u > 0) {
                        if ((last_add != u || last_remove != i) && (p == -1 || (rand() & 1))) p = i, q = u;
                    }
                    if (score_v > 0) {
                        if ((last_add != u || last_remove != i) && (p == -1 || (rand() & 1))) p = i, q = v;
                    }
                }
            break;
        }
    }
    if (p == -1) {
        auto it = UL.begin();
        while(!UL.empty()) {
            int u = it->u, v = it->v;
            if (covered[u][v] || it->lastuncover != lastuncover[u][v]) {
                it = UL.erase(it);
                continue;
            }
            for (int i = 0; i < n; ++i)
                if (C[i]) {
                    int score_u = dscore[i] + dscore[u] + w[i][u];
                    int score_v = dscore[i] + dscore[v] + w[i][v];
                    if (score_u > 0) {
                        if (p == -1 || (rand() & 1)) p = i, q = u;
                    }
                    if (score_v > 0) {
                        if (p == -1 || (rand() & 1)) p = i, q = v;
                    }
                }
            it = UL.erase(it);
            if (p != -1) break;
        }
    }
    if (p == -1) {
        last_add = last_remove = -1;
        return false;
    }
    last_add = q, last_remove = p;
    addtoC(q);
    removefromC(p);
    return true;
}

void greedy_update() {
    static bool tmp[N];
    static int f[N];
    memcpy(tmp, C, sizeof C);
    memset(f, 0, sizeof f);
    for (int i = 0; i < n; ++i)
        if (!tmp[i])
            for (int j = 0; j < n; ++j)
                if (G[i][j] && !tmp[j]) f[i]++;
    while(true) {
        int k, mx = 0;
        for (int i = 0; i < n; ++i)
            if (!tmp[i] && f[i] > mx) k = i, mx = f[i];
        if (mx == 0) break;
        tmp[k] = true;
        for (int i = 0; i < n; ++i)
            if (G[k][i] && !tmp[i]) f[i]--;
    }
    memcpy(_C, tmp, sizeof tmp);
}

void MC(int delta, int maxSteps) {
    Csize = 0, Lsize = m;
    srand(19260817);
    for (int i = 0; i < n; ++i) {
        int cnt = 0;
        for (int j = 0; j < n; ++j)
            if (i!=j && G[i][j]) {
                ++cnt, w[i][j]=1;
                if (i < j) {
                    L.push_back((edge){0,i,j});
                    UL.push_back((edge){0,i,j});
                }
            } 
        dscore[i] = cnt;
    }
    for (int i = 0; i < n; ++i) {
        int k = -1, mx = 0;
        for (int j = 0; j < n; ++j)
            if (!C[j] && dscore[j] >= mx) {
                if (dscore[j] == mx && (rand() & 1)) continue;
                mx = dscore[j], k = j;
            }
        if (mx == 0) break;
        addtoC(k);
    }
    int ub = Csize;
    memcpy(_C, C, sizeof (C));
    remove_some(ub, delta);

    int last_add = -1, last_remove = -1;
    for (step = 0; step < maxSteps; ++step) {
        if (!ChooseExchangePair(last_add, last_remove)) {
            for (auto it = L.begin(); it != L.end();) {
                int u = it->u, v = it->v;
                if (covered[u][v] || it->lastuncover != lastuncover[u][v])
                    it = L.erase(it);
                else {
                    w[u][v]++, w[v][u]++;
                    dscore[u]++, dscore[v]++;
                    ++it;
                }
            }
            int uindex = rand()%Csize, vindex = rand()%(n-Csize);
            int ucnt = 0, vcnt = 0;
            int u, v;
            for (int i = 0; i < n; ++i) {
                if (C[i] && (ucnt++) == uindex) u = i;
                if (!C[i] && (vcnt++) == vindex) v = i;
            }
            addtoC(v);
            removefromC(u);
        }

        if (Csize + Lsize < ub) {
            ub = Csize + Lsize;
            greedy_update();
        }
        remove_some(ub, delta);
    }

    int cnt = 0;
    for (int i = 0 ; i < n; ++i)
        if (!_C[i]) cnt++;
    printf("%d\n", cnt);
    int i;
    for (i = 0; i < n; ++i)
        if (!_C[i]) {
            printf("%d", i+1);
            break;
        }
    for (++i; i < n; ++i)
        if (!_C[i]) printf(" %d", i+1);
    puts("");
}

int main() {
    //freopen("in.txt","r",stdin);
    while (~scanf("%d%d",&n,&m)) {
    init();
    for (int i = 0,u,v; i < m; ++i) {
        scanf("%d%d",&u,&v); u--; v--;
        if (u != v) G[u][v] = G[v][u] = true;
    }
    m = 0;
    for (int i = 0; i < n; ++i)
        for (int j = i+1; j < n; ++j) {
            if (G[i][j]) G[i][j] = G[j][i] = 0;
            else {
                m++;
                G[i][j] = G[j][i] = 1;
            }
        }
    MC(2, 100000);
    }
    return 0;
}
