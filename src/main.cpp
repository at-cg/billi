#include<iostream>
#include<fstream>
#include<filesystem>
#include<vector>
#include<map>
#include<algorithm>
#include<stack>
#include<regex>
#include<cassert>
#define ll long long int
#define pb push_back
#define rb pop_back
#define ti tuple<int, int, int>
#define pii pair<int, int>
#define piii pair<int, pii>
#define piiii pair<int, piii>
#define mp make_pair
#define mt make_tuple
#define F first
#define S second

using namespace std;
namespace fs = std::filesystem;

fstream f;
string inputpath, outputdir;
int n, maxd; // max depth possible = total number of nodes
// vector<int> deg;

struct edge{
    int low, high; // a bracket can uniquely be identified by the lower and higher vertices it connects to : low -> lower height
    edge* front;
    edge* back; // for doubly linked list implementation -> helps in O(1) deletion
    edge(int l, int h):low(l), high(h), front(nullptr), back(nullptr){}
};

struct bracketlist{
    int sz, d1, d2; // d1 -> first minimum depth, d2 -> second minimum depth
    edge* start;
    edge* end; // will need the end as well for merging and finding the topmost bracket
    bracketlist():sz(0), d1(maxd), d2(maxd), start(nullptr), end(nullptr){}
    bracketlist(int lsz, int depth1, int depth2, edge* pstart, edge* pend): sz(lsz), d1(depth1), d2(depth2), start(pstart), end(pend){}
};

vector<vector<int>> g;
vector<vector<edge*>> remove_brackets;
vector<bool> mark; // for marking back edges
vector<pii> canonical_sese; // for storing the canonical sese pairs
vector<int> stack_trace; // for storing the vertices in the stack during dfs traversal
vector<int> depth; // for storing the depth of the vertices in the spanning tree
stack<piiii> st; // for storing the <top_bracket, list_size> for the black edges

map<string, int> lmap; // for storing the label for a particular gene
vector<string> ilmap; // for storing the gene for a particular label

int Ss = -1, Se = -1; // for connecting the tips

inline int last_bit(string s){
    return s == "+" ? 0 : 1;
}

inline string get_label(int x){
    return ilmap[x >> 1] + " " + (x & 1 ? "-" : "+");
}

void get_n(){
    f.open(inputpath, ios::in);
    string line;
    regex strip("^\\s+|\\s+$"), split("\\t");
    vector<string> tokens;
    n = 0; // number of nodes in the pangene graph
    if(f.is_open()){
        while(getline(f, line)){
            tokens.clear();
            line = regex_replace(line, strip, "");        
            sregex_token_iterator it(line.begin(), line.end(), split, -1), end;
            while(it != end){
                tokens.pb(*it);
                it++;
            }
            if(tokens[0] == "S"){
                assert(tokens.size() >= 3);
                lmap[tokens[1]] = n;
                ilmap.pb(tokens[1]);
                n++;
            }else if(tokens[0] == "L"){
                f.close();
                return;
            }
        }
        f.close();
    }
}

void make_graph(){
    f.open(inputpath, ios::in);
    string line;
    regex strip("^\\s+|\\s+$"), split("\\t");
    vector<string> tokens;
    if(f.is_open()){
        while(getline(f, line)){
            tokens.clear();
            line = regex_replace(line, strip, "");        
            sregex_token_iterator it(line.begin(), line.end(), split, -1), end;
            while(it != end){
                tokens.pb(*it);
                it++;
            }
            if(tokens[0] == "L"){
                assert(tokens.size() >= 5);
                int n1 = lmap[tokens[1]]; string s1 = tokens[2]; 
                int n2 = lmap[tokens[3]]; string s2 = tokens[4]; 
                // n1 and n2 are 0-indexed
                int id1 = (n1 << 1) + last_bit(s1), id2 = (n2 << 1) + last_bit(s2);
                // deg[(n1 << 1) + last_bit(s1)]++; deg[(n2 << 1) + last_bit(s2)]++;
                g[id1].pb(id2); g[id2].pb(id1);
            }
        }
        f.close();
    }
}

void merge(bracketlist* bl1, bracketlist* bl2){
    if(!bl2->start)return; // no benefit of merging
    bl1->end = bl2->end;
    if(!bl1->start){
        bl1->start = bl2->start;
        bl1->d1 = bl2->d1;
        bl1->d2 = bl2->d2;
    }else{
        bl1->end->front = bl2->start;
        bl2->start->back = bl1->end;
        if(bl1->d1 < bl2->d1){
            bl1->d2 = min(bl1->d2, bl2->d1);
        }else{
            bl1->d2 = min(bl1->d1, bl2->d2);
            bl1->d1 = bl2->d1;
        }
    }
    bl1->sz += bl2->sz;
    delete(bl2);
}

bracketlist* sese(int u, int parent){
    stack_trace.pb(u);
    mark[u] = true;
    if(parent != -1)depth[u] = depth[parent] + 1;

    bracketlist* bl = new bracketlist(); 
    for(int v : g[u]){
        if(v == parent)continue;
        bracketlist* bl1;
        if(mark[v]){// back-edge -> will come first in the dfs traversal for the node at greater depth
            if(depth[v] > depth[u])continue; // front-edge
            edge* ed = new edge(v, u);
            bl1 = new bracketlist(1, depth[v], maxd, ed, ed);
            remove_brackets[v].pb(ed);
        }else bl1 = sese(v, u);
        if((u ^ v) == 1){// check for canonical sese if it is a black edge, won't be true for back edge
            piiii val = st.top();
            if(bl1->sz > 0){
                piii cval = {bl1->sz, {bl1->end->low, bl1->end->high}};
                if(cval != get<1>(val)){
                    st.push({u, cval}); // u because it is the close end of the bibubble end
                }else{
                    canonical_sese.pb({v, get<0>(val)});// v because it is the close end of the bibubble start
                    st.pop();
                }
            }
        }
        merge(bl, bl1); // merging brackets from the child nodes
    }
    // removing brackets from respective lists after computing the value for the outgoing edges from that node
    for(edge* ed : remove_brackets[u]){
        if(ed->front){
            if(ed->back){
                ed->back->front = ed->front;
                ed->front->back = ed->back;
            }else{// is the first bracket
                ed->front->back = nullptr;
            }
        }else{// is the last bracket
            bl->end = ed->back;
            if(ed->back)ed->back->front = nullptr;
        }
        delete(ed);
        bl->sz--; // The current bracket list will only contain the edges that are being removed
    }
    // only left with the brackets that go up the node u
    // adding capping backedge
    // edge to be added u -> stack_trace[bl->d2]
    int w = stack_trace[bl->d2]; // stack_trace is 0-indexed and so are bl->d1 and bl->d2
    // edge need not be added to graph g
    // check if the edge already exists
    bool exists = false;
    for(int v : g[u]){
        if(v == w){
            exists = true; break;
        }
    }
    if(!exists){
        edge* ed = new edge(w, u);
        bracketlist* bl1 = new bracketlist(1, depth[w], maxd, ed, ed);
        remove_brackets[w].pb(ed);
        merge(bl, bl1);
    }

    stack_trace.rb();
    return bl;
}

int main(int argc, char* argv[])
{   
    ios_base::sync_with_stdio(false); cin.tie(0); cout.tie(0); // Fast IO
    if(argc != 3){
        return 1;
    }
    inputpath = argv[1];
    outputdir = argv[2];
    if(!fs::exists(outputdir))fs::create_directories(outputdir);
    string summarypath = outputdir + "/summary.txt";
    freopen(summarypath.c_str(), "w", stdout);

    get_n();
    // deg.resize(2 * n);

    // for a node x (0-indexed) in pangene graph - two nodes 2 * x (+), 2 * x + 1 (-) are created in the bi-edged graph 
    // using vectors will be good coz we will have to add capping backedges as well
    g.resize(2 * n); // // +2 is for S if required -> not needed
    make_graph(); // adding grey edges
    // adding black edges
    for(int i = 0; i < n; i++){
        int n1 = 2 * i, n2 = 2 * i + 1;
        g[n1].pb(n2); g[n2].pb(n1);
    }

    // Graph compression - removing unitigs --- Not implementing currently -> does not seem that useful

    // Identifying tips
    // Basically the vertices with deg = 1
    vector<int> tips; // vector for storing the vertices that represent the tip
    for(int i = 0; i < 2 * n; i++){
        if(g[i].size() == 1){
            tips.pb(i);
            assert((g[i][0] ^ i) == 1); // checking whether the node is connected via a black edge
        }
    }

    if(tips.size() == 0){ // there can be multiple components here, need to take care of in the dfs 
        Ss = 0; Se = 1; // no extra edges required
    }else{
        Ss = tips[0]; Se = Ss ^ 1; 
        // if(tips.size() > 1){// else no need to add an extra edge
            for(int i = 1; i < tips.size(); i++){
                g[Ss].pb(tips[i]); g[tips[i]].pb(Ss);
            }
        // }
    }

    // SESE
    // An important observation is that bibubble ends won't turn as backedges
    mark.resize(2 * n);
    remove_brackets.resize(2 * n);
    depth.resize(2 * n);
    
    maxd = 2 * n + 100;
    sese(Ss, -1); 
    assert(st.size() == 0);

    // Taking care of multiple components
    for(int i = 0; i < 2 * n; i += 2){
        if(mark[i]){
            assert(mark[i + 1]);
            continue;
        }
        Ss = i; Se = i + 1;
        sese(Ss, -1);
        assert(st.size() == 0);
    }

    for(pii rs : canonical_sese)cout << get_label(rs.F) << " " << get_label(rs.S) << endl;
} 
