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
    int id; // a bracket can uniquely be identified by the lower and higher vertices it connects to : low -> lower height, assigning a unique to each
    edge* front;
    edge* back; // for doubly linked list implementation -> helps in O(1) deletion
    edge(int eid):id(eid), front(nullptr), back(nullptr){}
};

struct bracketlist{
    int sz, d1, d2; // d1 -> first minimum depth, d2 -> second minimum depth
    edge* start;
    edge* end; // will need the end as well for merging and finding the topmost bracket
    bracketlist():sz(0), d1(maxd), d2(maxd), start(nullptr), end(nullptr){}
    bracketlist(int lsz, int depth1, int depth2, edge* pstart, edge* pend): sz(lsz), d1(depth1), d2(depth2), start(pstart), end(pend){}
};

template <typename... Args>
void printArgs(Args... args){
    return;
    ((cout << args << " "), ...) << endl; // Fold expression
}

void printBracketList(bracketlist* bl){
    return;
    printArgs("d1:", bl->d1, "d2:", bl->d2);
    cout << "Bracket List: ";
    edge* it = bl->start;
    while(it){
        cout << it->id << " ";
        it = it->front;
    }cout << endl;
}

// Not required : last_bit, upd variable in merge function

vector<vector<pii>> g;
vector<vector<edge*>> remove_brackets;
vector<bool> mark; // for marking back edges
vector<pii> canonical_sese; // for storing the canonical sese pairs
vector<int> stack_trace; // for storing the vertices in the stack during dfs traversal
vector<int> depth; // for storing the depth of the vertices in the spanning tree
vector<pii> eid; 
vector<pii> st; // (node, size) will store the size of the bracket list when that bracket is the topmost bracket and for which node that was topmost previously
map<string, int> lmap; // for storing the label for a particular gene
vector<string> ilmap; // for storing the gene for a particular label

int Ss = -1, Se = -1; // for connecting the tips
int tot_grey = 0; // for storing total number of grey edges, will help in accessing the brackets

// inline int last_bit(string s){
//     return s == "+" ? 0 : 1;
// }

inline string get_label(int x, int y){
    return (ilmap[x >> 1] + " " + (x & 1 ? "+" : "-") + " " + ilmap[y >> 1] + " " + (y & 1 ? "-" : "+"));
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
                int id1 = n1 << 1, id2 = n2 << 1;
                if(s1 == "+")id1++;
                if(s2 == "-")id2++;
                // deg[(n1 << 1) + last_bit(s1)]++; deg[(n2 << 1) + last_bit(s2)]++;
                if(n1 == n2)continue; // there will be a black edge b/w them <- safe operation
                tot_grey++;
                // printArgs(id1, id2, tot_grey);
                g[id1].pb({id2, tot_grey}); g[id2].pb({id1, tot_grey});
                eid.pb({min(id1, id2), max(id1, id2)});
            }
        }
        f.close();
    }
}

void merge(bracketlist* bl1, bracketlist* bl2, bool upd){
    if(!bl1->start){
        bl1->start = bl2->start;
        if(upd){
            bl1->d1 = bl2->d1;
            bl1->d2 = bl2->d2;
        }
    }else{
        bl1->end->front = bl2->start;
        bl2->start->back = bl1->end;
        if(upd){
            if(bl1->d1 < bl2->d1){
                bl1->d2 = min(bl1->d2, bl2->d1);
            }else{
                bl1->d2 = min(bl1->d1, bl2->d2);
                bl1->d1 = bl2->d1;
            }
        }
    }
    bl1->end = bl2->end;
    bl1->sz += bl2->sz;
    delete(bl2);
}

bracketlist* sese(int u, int parent){
    stack_trace.pb(u);
    mark[u] = true;
    if(parent != -1)depth[u] = depth[parent] + 1;

    printArgs(u, parent); 
    bracketlist* bl = new bracketlist(); 
    int cnt_back = 0;
    for(pii child : g[u]){
        int v = child.F;
        if(v == parent)continue;
        bracketlist* bl1;
        if(mark[v]){// back-edge -> will come first in the dfs traversal for the node at greater depth
            if(depth[v] > depth[u])continue; // front-edge
            assert(child.S != -1);
            edge* ed = new edge(child.S - 1); // tot_grey count starts from 1
            bl1 = new bracketlist(1, depth[v], maxd, ed, ed);
            remove_brackets[v].pb(ed);
        }else bl1 = sese(v, u);
        if((u ^ v) == 1){// check for canonical sese if it is a black edge, won't be true for back edge
            printArgs("Bracket list size on reaching the solid edge:", bl1->sz);
            if(bl1->sz > 0){
                int br = bl1->end->id;
                assert(br != -1);
                if(st[br] != mp(0, 0)){
                    printArgs("Found the canonical pair");
                    assert(st[br].S == bl1->sz);
                    if(g[v].size() > 2 || g[st[br].F].size() > 2){ // to avoid linear chains
                        printArgs("Found the contributing canonical pair:", v, st[br].F);
                        canonical_sese.pb({v, st[br].F});
                    }
                }
                st[br] = {u, bl1->sz}; // will happen regardless you found something or not
            }
        }
        if(bl1->start){
            printArgs(u, "child", v);
            printBracketList(bl); printBracketList(bl1);
            merge(bl, bl1, true); // merging brackets from the child nodes
            printBracketList(bl);
            if(bl1->d1 < depth[u] && depth[v] > depth[u])cnt_back++; // an edge from u to its ancestor won't contribute to capping backedge
        }
    }

    // removing brackets from respective lists after computing the value for the outgoing edges from that node
    for(edge* ed : remove_brackets[u]){
        printArgs("remove", eid[ed->id].F, eid[ed->id].S);
        int min_depth = min(depth[eid[ed->id].F], depth[eid[ed->id].S]);
        if(bl->d1 == min_depth)bl->d1 = maxd;
        else if(bl->d2 == min_depth)bl->d2 = maxd;
        if(ed->front){
            if(ed->back){
                ed->back->front = ed->front;
                ed->front->back = ed->back;
            }else{// is the first bracket
                bl->start = ed->front;
                ed->front->back = nullptr;
            }
        }else{// is the last bracket
            bl->end = ed->back;
            if(ed->back)ed->back->front = nullptr;
            else{
                bl->start = bl->end = nullptr;
            }
        }
        delete(ed);
        bl->sz--; // The current bracket list will only contain the edges that are being removed
        printBracketList(bl);
    }
   
    // only left with the brackets that go up the node u
    // adding capping backedge
    // edge to be added u -> stack_trace[bl->d2]
    if(cnt_back >= 2){// at least two childs should have backedges with depth lower than node u in order to add a capping backedge
        int w = stack_trace[bl->d2]; // stack_trace is 0-indexed and so are bl->d1 and bl->d2
        // edge need not be added to graph g
        // check if the edge already exists
        bool exists = false;
        for(pii child : g[u]){
            if(child.F == w){
                exists = true; break;
            }
        }
        printArgs("Capping back-edge", u, w);
        if(!exists){
            printArgs("Capping back-edge does not exist");
            tot_grey++;
            edge* ed = new edge(tot_grey);
            eid.pb({min(w, u), max(w, u)});
            bracketlist* bl1 = new bracketlist(1, depth[w], maxd, ed, ed);
            remove_brackets[w].pb(ed);
            merge(bl, bl1, true);
        }
    }

    // printArgs(u, parent, stack_trace.size(), bl->sz);
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
    // adding black edges first <- important since don't want black edge to appear as a back or front edge
    for(int i = 0; i < n; i++){
        int n1 = 2 * i, n2 = 2 * i + 1;
        g[n1].pb({n2, -1}); g[n2].pb({n1, -1});
    }
    make_graph(); // adding grey edges
    
    // Graph compression - removing unitigs --- Not implementing currently -> does not seem that useful

    // Identifying tips
    // Basically the vertices with deg = 1
    vector<int> tips; // vector for storing the vertices that represent the tip
    for(int i = 0; i < 2 * n; i++){
        if(g[i].size() == 1){
            tips.pb(i);
            assert((g[i][0].F ^ i) == 1); // checking whether the node is connected via a black edge
        }
    }

    // cout << "TIPS\n";
    // for(int x : tips)cout << x << " ";
    // cout << endl;

    if(tips.size() == 0){ // there can be multiple components here, need to take care of in the dfs 
        Ss = 0; Se = 1; // no extra edges required
    }else{
        Ss = tips[0]; Se = Ss ^ 1; 
        // if(tips.size() > 1){// else no need to add an extra edge
            for(int i = 1; i < tips.size(); i++){
                tot_grey++;
                g[Ss].pb({tips[i], tot_grey}); g[tips[i]].pb({Ss, tot_grey}); // edge id's are required here as well
                eid.pb({min(Ss, tips[i]), max(Ss, tips[i])});
            }
        // }
    }

    // for(int i = 0; i < eid.size(); i++){
    //     printArgs(i, eid[i].F, eid[i].S);
    // }
    
    // SESE
    // An important observation is that bibubble ends won't turn as backedges
    mark.resize(2 * n);
    remove_brackets.resize(2 * n);
    depth.resize(2 * n);
    st.resize(tot_grey + n + 1); // grey_edges count start from 1 and every vertex will have at most one back edge

    maxd = 2 * n + 100;
    sese(Ss, -1);     
    assert(stack_trace.size() == 0);
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
    cout << "Total canonical cycle equivalent pairs found: " << canonical_sese.size() << endl;
    for(pii rs : canonical_sese)cout << get_label(rs.F, rs.S) << endl;
    f.flush();
    f.close();
} 
