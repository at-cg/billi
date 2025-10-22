#include<CLI11.hpp>
#include<iostream>
#include<fstream>
#include<filesystem>
#include<vector>
#include<map>
#include<set>
#include<unordered_map>
#include<algorithm>
#include<stack>
#include<queue>
#include<deque>
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
#define PRINT true

using namespace std;
namespace fs = std::filesystem;

// ****************************************************************************************************************************************************** 
// **** Generic ****
// ****************************************************************************************************************************************************** 
const ll maxv = __LONG_LONG_MAX__;
int tot_grey = 0; // for storing total number of grey edges, will help in accessing the brackets
int maxd; // max depth possible = total number of nodes
string summarypath; // path where the given stats will be written
vector<bool> mark; // for marking the vertices that are have been visited
// ****************************************************************************************************************************************************** 

// ****************************************************************************************************************************************************** 
// **** Structures ****
// ****************************************************************************************************************************************************** 
struct edge{
    int id; // a bracket can uniquely be identified by the lower and higher vertices it connects to : low -> lower height, assigning a unique to each
    int node; // required node corresponding to the black edge when edge 'id' is the top in the bracketlist
    int sz; // required size of bracketlist corresponding to the black edge when edge 'id' is the top in the bracketlist
    edge* front;
    edge* back; // for doubly linked list implementation -> helps in O(1) deletion
    edge(int eid):id(eid), node(-1), sz(-1), front(nullptr), back(nullptr){}
};

struct bracketlist{
    int sz, d; // d -> minimum depth
    edge* start;
    edge* end; // will need the end as well for merging and finding the topmost bracket
    bool merged = false;
    bracketlist(): sz(0), d(maxd), start(nullptr), end(nullptr){}
    bracketlist(int lsz, int depth, edge* pstart, edge* pend): sz(lsz), d(depth), start(pstart), end(pend){}
};
// ****************************************************************************************************************************************************** 

// ****************************************************************************************************************************************************** 
// **** Input ****
// ****************************************************************************************************************************************************** 
fstream f;
string inputpath, outputdir;
bool print_allele = false, print_equivalent = false, print_hairpin = false, print_panbubble_tree = false; // whether to report alleles, cycle equivalent pairs and hairpins
int maxsize = 1000, maxdepth = 1; // max size of panbubbles whose vertices are also to be printed, max depth to which the panbubbles are to be reported
int n, edges; // no of nodes (genes), no of edges
// int before_initialisation_comp_cnt = 0; // no of connected components in the input graph before initialisation
vector<string> ilmap; // for storing the gene for a particular label
vector<bool> has_self_loop; // whether the node has a self loop <- can't be cycle equivalent with any other node
vector<vector<pii>> g; // (node_id, grey_edge_id)
map<string, int> lmap; // for storing the label for a particular gene
map<int, int> id_in_original_graph; // used to find the id in original graph, in order to use the property of black edges map[u] ^ map[v] = 1

void get_ne() {
    ifstream f(inputpath);
    string line;
    n = 0;
    edges = 0;

    if (f.is_open()) {
        while (getline(f, line)) {
            // Trim leading and trailing whitespace manually
            size_t start = line.find_first_not_of(" \t\r\n");
            size_t end = line.find_last_not_of(" \t\r\n");

            if (start == string::npos || end == string::npos)
                continue; // Skip empty or all-whitespace lines

            line = line.substr(start, end - start + 1);

            // Tokenize manually using '\t' as delimiter
            vector<string> tokens;
            size_t pos = 0, prev = 0;
            while ((pos = line.find('\t', prev)) != string::npos) {
                tokens.emplace_back(line.substr(prev, pos - prev));
                prev = pos + 1;
            }
            tokens.emplace_back(line.substr(prev)); // Last token

            if (!tokens.empty()) {
                if (tokens[0] == "S") {
                    lmap[tokens[1]] = n; // 0-indexed
                    ilmap.push_back(tokens[1]);
                    ++n;
                } else if (tokens[0] == "L") {
                    ++edges;
                } //else if (n != 0) break;
            }
        }
        f.close();
    }
}

void add_edge(vector<vector<pii>>& g, int id1, int id2){
    g[id1].pb({id2, tot_grey}); 
    if(id1 != id2) g[id2].pb({id1, tot_grey}); // 0-indexed
    tot_grey++;
}

void make_graph() {
    ifstream f(inputpath);
    string line;

    if (f.is_open()) {
        while (getline(f, line)) {
            // Manual trim (faster than regex)
            size_t start = line.find_first_not_of(" \t\r\n");
            size_t end = line.find_last_not_of(" \t\r\n");
            if (start == string::npos || end == string::npos)
                continue;
            line = line.substr(start, end - start + 1);

            // Manual tab-based tokenization
            vector<string> tokens;
            size_t prev = 0, pos;
            while ((pos = line.find('\t', prev)) != string::npos) {
                tokens.emplace_back(line.substr(prev, pos - prev));
                prev = pos + 1;
            }
            tokens.emplace_back(line.substr(prev));  // Last token

            // Process only 'L' lines
            if (!tokens.empty() && tokens[0] == "L") {
                int n1 = lmap[tokens[1]];
                string s1 = tokens[2];
                int n2 = lmap[tokens[3]];
                string s2 = tokens[4];

                // Consistency of edge direction
                if (n1 > n2) {
                    swap(n1, n2);
                    s1 = (s1 == "+") ? "-" : "+";
                    s2 = (s2 == "+") ? "-" : "+";
                    swap(s1, s2);
                }

                int id1 = n1 << 1;
                int id2 = n2 << 1;
                if (s1 == "+") id1++;
                if (s2 == "-") id2++;

                add_edge(g, id1, id2);

                if (n1 == n2 && s1 != s2) {
                    has_self_loop[id1] = true;
                }
            }
        }
        f.close();
    }
}

// ****************************************************************************************************************************************************** 

// ****************************************************************************************************************************************************** 
// **** Helper functions ****
// ****************************************************************************************************************************************************** 
string get_single_label(int x, int ty){
    x = id_in_original_graph[x];
    return ((ty == 0 ? (x & 1 ? ">" : "<") : (x & 1 ? "<" : ">")) + ilmap[x >> 1]);
}

template <typename... Args>
void printArgs(Args... args){
    if(!PRINT)return;
    ((cout << args << " "), ...) << endl; // Fold expression
}

void printBracketList(bracketlist* bl){
    if(!PRINT)return;
    // printArgs("d1:", bl->d1, "d2:", bl->d2);
    cout << "Bracket List: ";
    edge* it = bl->start;
    while(it){
        cout << it->id << " ";
        it = it->front;
    }cout << endl;
}

void printVector(vector<int>& v){
    if(!PRINT)return;
    for(int x : v){
        cout << x << " ";
    }
    cout << endl;
}

void printVector(vector<pii>& v, int ty = 0){
    if(!PRINT)return;
    for(pii x : v){
        cout << "(" << (ty == 0 ? to_string(x.F) : get_single_label(x.F, 0)) << ", " << x.S << ") ";
    }
    cout << endl;
}

void printVector(vector<char>& v){
    if(!PRINT)return;
    for(char x : v){
        cout << x << " ";
    }
    cout << endl;
}

void printGraph(vector<vector<int>>& g){
    if(!PRINT)return;
    int sz = g.size();
    for(int i = 0; i < sz; i++){
        cout << i << ": ";
        printVector(g[i]);
    }
}

void printGraph(vector<vector<pii>>& g){
    if(!PRINT)return;
    int sz = g.size();
    for(int i = 0; i < sz; i++){
        cout << get_single_label(i, 0) << ": ";
        printVector(g[i], 1);
    }
}
// ****************************************************************************************************************************************************** 

// ****************************************************************************************************************************************************** 
// **** Processed Graph ****
// ****************************************************************************************************************************************************** 
int after_initialisation_comp_cnt = 0; // no of connected components in the input graph after initialisation
int id_ptr; // ptr to the id number from which new nodes will be added
int n_processed; // no of nodes in the processed graph (at most {n + 2 * edges})
vector<char> type_edge; // finding the type of edge split at the black edge that is a bridge (0 -> no split, 1 -> split from u, 2 -> split from v, 3 -> split from both, 4 -> edge not of any use)
vector<int> first_node; // for storing a node in every connected-component in the biedged graph
vector<int> biedged_connected_comp; // connected-component in the biedged graph to which a node belongs to 
vector<int> backedge_cnt; // upper bound on number of back edges
vector<int> Start; // starting point of search for each component in original graph
vector<int> depth; // for storing the depth of the vertices in the spanning tree
vector<int> dual; // id of the vertex connected via a black edge to the given vertex
vector<ll> tip_start; // will store the id from which extra edges because of tips are added for each biedged_connected_comp
vector<vector<int>> tips; // vector for storing the vertices that represent the tip
vector<int> tot_grey_comp; // count of grey edges in the processed graph for individual components
vector<vector<pii>> g_processed; // (node_id, grey_edge_id)

set<int> find_unique_excluding_node(int u, int v){// exclude node u and v
    set<int> s;
    for(pii child : g[u]){
        if(child.F == u || child.F == v)continue;
        s.insert(child.F);
        if(s.size() == 2)break;
    }
    return s;
}

void neigh_copy(int u, int x){
    // black edges added in the start
    if(x != -1){
        int v = id_ptr++; // new node created
        id_in_original_graph[v] = x;
        g_processed[v].pb({u, -1}); // will be a black edge, tot_grey won't be incremented
        g_processed[u].pb({v, -1}); 
        // edges += 2; // not to count the black edges
        dual[u] = v; dual[v] = u;
    }else{
        if(u & 1){// done only once for the odd labeled vertex
            dual[u] = u ^ 1; dual[u ^ 1] = u;
        }
    }
    // no new grey edges added
    for(pii child : g[u]){
        if(child.F == x)continue;
        g_processed[u].pb(child);
        // edges++;
    }
    id_in_original_graph[u] = u;
}

void dfs_comp(int u){ // Finding connected components in the biedged graph
    mark[u] = true;
    biedged_connected_comp[u] = after_initialisation_comp_cnt;
    for(pii child : g_processed[u]){
        int v = child.F;
        if(mark[v])continue;
        dfs_comp(v);
    }
}

int find_unique_including_selfloop(int u){ // skip u ^ 1
    set<int> s;
    for(pii child : g_processed[u]){
        if((id_in_original_graph[child.F] ^ id_in_original_graph[u]) == 1)continue;
        s.insert(child.F);
        break;
    }
    return s.size();
}
// ****************************************************************************************************************************************************** 

// ******************************************************************************************************************************************************  
// **** Finding bridges ****
// ******************************************************************************************************************************************************  
int timer = 0; // for finding bridges in the uninitialised graph
vector<int> tin, low;
vector<bool> mark_bridge;

void dfs_bridge(int u, int parent){
    // printArgs("dfs_bridge:", u, parent);
    mark[u] = true;
    tin[u] = low[u] = timer++;

    for(pii child : g[u]){
        int v = child.F;
        if(v == parent)continue; // a multiedge b/w 2 nodes can also result in a bridge
        if(mark[v])low[u] = min(low[u], tin[v]);
        else{
            dfs_bridge(v, u);
            low[u] = min(low[u], low[v]);
            // printArgs(u, v, low[v], tin[u]);
            if(low[v] > tin[u] && ((u ^ v) == 1)){
                // printArgs("Bridge found:", u, v);
                mark_bridge[u >> 1] = true;
            } // strictly greater sign (>) -- loop case
        }
    }
}
// ****************************************************************************************************************************************************** 

// ****************************************************************************************************************************************************** 
// **** Single Entry Single Exit Algorithm (Johnson et al., 1994) ****
// ******************************************************************************************************************************************************  
int nodes = 0; // no of nodes in a given component in the input graph
// int st0 = 0; // top node when the bracketlist is empty
// int time = 0;
ll multiplier = 1; // for computing the hash
int possible_pairs = 0, valid_pairs = 0, hairpins = 0; // storing the number of valid panbubble pairs

vector<int> stack_trace; // for storing the vertices in the stack during dfs traversal
vector<int> bubble_depth; // depth of the bubbles in the bubble tree
vector<int> rm_cnt; // vector storing the number of brackets ending at that node
// vector<int> bracket_seq; // for storing the id of the bracket sequences so as to make the bubble tree
// vector<int> visit_time; // for storing the visit time of the vertices during the dfs traversal
vector<bool> root_bg; // root bubble id's in the bubble tree
vector<bool> inb; // whether the node has been added in the bubble tree
vector<pii> canonical_sese; // for storing the canonical sese pairs
vector<vector<int>> bg; // vector for storing the adjacency matrix for the bubble graph
// unordered_map<ll, int> st; // (edge, size) -> hashed into a ll key ** Takes a lot of time **
map<ll, int> st; // (edge, size) -> hashed into a ll key
vector<bracketlist*> bl; // vector storing the bracket list for different nodes
vector<vector<edge*>> remove_brackets; // vector that stores the brackets to be removed while exiting a particular vertex during SESE algorithm
deque<int> qbg; // queue maintained while doing bfs on the bubble tree

ll get_key(int id, int size){
    return id * multiplier + size;
}

int find_unique_excluding_selfloop(int u, int v){// exclude nodes u, u ^ 1, v
    set<int> s;
    for(pii child : g_processed[u]){
        if((id_in_original_graph[child.F] ^ id_in_original_graph[u]) <= 1 || child.F == v)continue;
        s.insert(child.F);
        break;
    }
    return s.size();
}

void merge(bracketlist* bl1, bracketlist* bl2){
    if(!bl2->start)return;

    if(!bl1->start){
        bl1->start = bl2->start;
        bl1->d = bl2->d;
    }else{
        bl1->end->front = bl2->start;
        bl2->start->back = bl1->end;
        bl1->d = min(bl1->d, bl2->d);
        // if(bl1->d1 < bl2->d1){
        //     bl1->d2 = min(bl1->d2, bl2->d1);
        // }else if(bl1->d1 > bl2->d1){
        //     bl1->d2 = min(bl1->d1, bl2->d2);
        //     bl1->d1 = bl2->d1;
        // }else{
        //     bl1->d2 = min(bl1->d2, bl2->d2);
        // }
        // else{
        //     bl1->d2 = min(bl1->d1, bl2->d2);
        //     bl1->d1 = bl2->d1;
        // }
    }
    bl1->end = bl2->end;
    bl1->sz += bl2->sz;
    // delete(bl2);
    bl2->merged = true;
}

void dfs_depth(int u, int parent, int comp){
    // assert(!mark[u]);
    if(parent != -1)depth[u] = depth[parent] + 1;
    nodes++;
    mark[u] = true;
    for(pii child : g_processed[u]){
        int v = child.F;
        if(mark[v]){
            if(depth[v] < depth[u] && v != parent)backedge_cnt[comp]++;
            continue;
        }
        dfs_depth(v, u, comp);
    }
}

void reinitialise_edge_ids(int u, int parent, int& comp, int& pass){
    mark[u] = true;
    for(pii child : g_processed[u]){
        int v = child.F;
        if(mark[v])continue;
        if(child.S != -1){
            if(pass == 1 || (pass == 0 && child.S < tip_start[comp])){
                if(pass == 1 && child.S == tip_start[comp]){
                    tip_start[comp] = child.S = tot_grey_comp[comp]++;
                }else child.S = tot_grey_comp[comp]++;
            }
        }
        reinitialise_edge_ids(v, u, comp, pass);
    }
}

void sese_minbracket(int u, int parent, int& x, int& val){
    // printArgs(u, parent, x, val);
    // assert(!mark[u]);
    mark[u] = true;
    if(parent != -1)depth[u] = depth[parent] + 1;

    for(pii child : g_processed[u]){
        int v = child.F;
        if(v == parent || v == u)continue; // v = u won't contribute anything and v == parent won't contribute to backedge
        if(mark[v]){// back-edge -> will come first in the dfs traversal for the node at greater depth
            if(depth[v] > depth[u])continue; // front-edge -- multi edges will be taken care of here
            // assert(child.S != -1); // black edges can't be backedges
        }else sese_minbracket(v, u, x, val);
    }

    bl[u] = new bracketlist(); 

    for(pii child : g_processed[u]){
        int v = child.F;
        if(depth[v] != depth[u] + 1)continue;
        if(bl[v] && !bl[v] -> merged){
            // printArgs(u, v, bl[v]->sz);
            if((id_in_original_graph[u] ^ id_in_original_graph[v]) == 1){
                int bridge_cnt = bl[v]->sz;
                if(bridge_cnt < val){
                    val = bridge_cnt;
                    x = u;
                }
            }

            bl[u]->sz += bl[v]->sz;
            bl[v]->merged = true;
        }
    }

    // removing brackets from respective lists
    bl[u]->sz -= rm_cnt[u];
    
    // pushing back edges from node u
    for(pii child : g_processed[u]){
        int v = child.F;
        if(depth[v] >= depth[u] - 1)continue;
        rm_cnt[v]++;
        bl[u]->sz++;
    }
}

void sese(int u, int parent){
    stack_trace.pb(u);
    mark[u] = true;

    int h0 = maxd, h1 = maxd, h2 = maxd;

    for(pii child : g_processed[u]){
        int v = child.F;
        if(v == parent || v == u)continue; // v = u won't contribute anything and v == parent won't contribute to backedge
        if(mark[v]){// back-edge -> will come first in the dfs traversal for the node at greater depth
            if(depth[v] > depth[u])continue; // front-edge -- multi edges will be taken care of here
            // assert(child.S != -1); // black edges can't be backedges
            h0 = min(h0, depth[v]);
        }else sese(v, u);
    }

    bl[u] = new bracketlist(); 

    for(pii child : g_processed[u]){
        int v = child.F;
        if(depth[v] != depth[u] + 1)continue;
        if(bl[v]){
            if(bl[v]->d < h1){
                h2 = h1;
                h1 = bl[v]->d;
            }else{
                h2 = min(h2, bl[v]->d);
            }
            bl[v]->d = maxd; // not to be used again
        }
        if(!bl[v] -> merged){
            // printBracketList(bl[u]); printBracketList(bl[v]);
            merge(bl[u], bl[v]);
            // printBracketList(bl[u]);
        }
    }
    bl[u]->d = min(h0, h1);

    // printArgs(h0, h1, h2);

    // removing brackets from respective lists
    for(edge* ed : remove_brackets[u]){
        if(ed->front){
            if(ed->back){
                ed->back->front = ed->front;
                ed->front->back = ed->back;
            }else{// is the first bracket
                bl[u]->start = ed->front;
                ed->front->back = nullptr;
            }
        }else{// is the last bracket
            bl[u]->end = ed->back;
            if(ed->back)ed->back->front = nullptr;
            else{
                bl[u]->start = nullptr;
                // bl->end = nullptr; // not required
            }
        }
        delete(ed);
        bl[u]->sz--; // The current bracket list will only contain the edges that are being removed
    }
    remove_brackets[u].clear();
    if(bl[u]->d == depth[u])bl[u]->d = maxd;

    // pushing back edges from node u
    for(pii child : g_processed[u]){
        int v = child.F;
        if(depth[v] >= depth[u] - 1)continue;
        edge* ed = new edge(child.S); // tot_grey count starts from 0
        remove_brackets[v].pb(ed);
        merge(bl[u], new bracketlist(1, depth[v], ed, ed));
    }

    // capping back edge
    if(h2 < h0){
        // printArgs("finding capping backedge");
        int w = stack_trace[h2];
        edge* ed = new edge(tot_grey); tot_grey++; // need not modify g
        remove_brackets[w].pb(ed);
        merge(bl[u], new bracketlist(1, depth[w], ed, ed));
        // printArgs("capping back edge:", u, w);
    }

    // finding cycle equivalence
    if(parent != -1 && ((id_in_original_graph[parent] ^ id_in_original_graph[u]) == 1)){
        // printArgs("finding cycle equivalence:", u, parent);
        ll key;
        if(bl[u]->sz > 0){
            int br = bl[u]->end->id;
            // assert(br != -1); // not a black edge
            key = get_key(br, bl[u]->sz);
        }else{
            key = 0;
        }

        if(st.find(key) != st.end()){
            int w = st[key];
            if(find_unique_excluding_selfloop(u, w) == 1 && find_unique_excluding_selfloop(w, u) == 1){
                // printArgs("cycle equivalent pair:", u, w);
                canonical_sese.pb({u, w});
            }
        }
        st[key] = parent; // will happen regardless you found something or not
    }

    stack_trace.rb();
}

// ******************************************************************************************************************************************************  

// ****************************************************************************************************************************************************** 
// **** SCC **** 
// ref - https://cp-algorithms.com/graph/strongly-connected-components.html extended to multigraph
// ******************************************************************************************************************************************************  
vector<bool> considered; // for saving which vertices have been considered twice while computing bb_comp
vector<bool> valid; // for marking whether a potential bi-bubble is valid
vector<int> bb_comp; // panbubble-component to which a node belongs to
vector<vector<int>> bb_nodes; // storing the nodes for different panbubbles
vector<int> opp_entrance; // for saving the other end of the panbubble
vector<int> aux_cc_comp; // connected-component in the auxiliary graph to which a node belongs to
vector<int> order; // will be a sorted list of G's vertices by exit time
vector<vector<int>> ag; // directed graph
vector<vector<int>> ag_rev; // create adjacency list of G^T   

string get_label(int x, int y){
    return (get_single_label(x, 0) + " " + get_single_label(y, 1));
}

string get_label_from_id(int id){
    return "(" + get_label(canonical_sese[id].F, canonical_sese[id].S) + ")";
}

bool end_gene(int& u, pii& rs){
    return u == rs.F || u == dual[rs.F] || u == rs.S || u == dual[rs.S];
}

void add_directed_edges_from_gene_endpoints(int id1, int id2){
    ag[dual[id1]].pb(id2); ag[dual[id2]].pb(id1);
}

void upd_node(int u, int id){
    mark[u] = true;
    if(bb_comp[u] != id) bb_nodes[id].pb(u);
    if(bb_comp[u] == -1)bb_comp[u] = id;
}

void mark_bb_nodes(int u, pii& rs, int id){
    int w = opp_entrance[u];
                
    if(w != -1){
        int v = dual[w];
        considered[w] = true;

        if(!end_gene(v, rs)){
            mark_bb_nodes(v, rs, id);
        }
        return;
    }

    // if((u == rs.S || u == dual[rs.S]) && bb_nodes[id].size() != 0)return; // u == dual[rs.S] -> can go to y complement but then the pair is not cycle equivalent, if size = 0 implies hairpin loop
    if(u != rs.F){// only the directed genes \in \tilde{U} set
        upd_node(u, id);
    }

    for(pii child : g_processed[u]){
        int v = child.F;
     
        if(child.S >= tip_start[biedged_connected_comp[rs.S]] || end_gene(v, rs))continue; // don't consider the edges added because of tips
        if(mark[v]){
            if(!considered[v]){    
                w = opp_entrance[v];
                if(w != -1){
                    v = dual[w];
                    considered[w] = true;
                }
                // if(bb_comp[v] != bb_comp[rs.F])bb_nodes[id].pb(v); // if the two panbubbles are adjacent, then too start and end points of the panbubble will be added, however the complete set need not be added                          
                if(!end_gene(v, rs)){
                    considered[v] = true;
                    mark_bb_nodes(v, rs, id);
                }
            }
            continue;
        }
        mark_bb_nodes(v, rs, id);
    }
}

void make_auxillary_graph_from_biedged_graph(){
    for(int i = 0; i < n_processed; i++){
        for(pii child : g_processed[i]){
            if(child.S >= tip_start[biedged_connected_comp[child.F]] || child.S == -1)continue; // don't consider the edges added because of tips and the black edges
            if(i <= child.F)add_directed_edges_from_gene_endpoints(i, child.F); // making sure edges are added only in one direction
        }
    }

    // using the extra added edges
    for(pii rs : canonical_sese){
        // printArgs("Overlapping edge:", dual[rs.F], dual[rs.S]);
        add_directed_edges_from_gene_endpoints(dual[rs.F], dual[rs.S]);
    }
}

// runs depth first search starting at vertex v.
// each visited vertex is appended to the output vector when dfs leaves it.
void scc_dfs(int v, vector<vector<int>> const& adj, vector<int>& output) {
    mark[v] = true;
    for (auto u : adj[v])
        if (!mark[u])
            scc_dfs(u, adj, output);
    output.push_back(v);
}

void scc(){
    fill(mark.begin(), mark.end(), false);

    // first series of depth first searches
    for (int i = 0; i < n_processed; i++)
        if (!mark[i])
            scc_dfs(i, ag, order);

    ag_rev.resize(n_processed);
    
    for (int v = 0; v < n_processed; v++)
        for (int u : ag[v])
            ag_rev[u].push_back(v);

    fill(mark.begin(), mark.end(), false);
    reverse(order.begin(), order.end());

    // second series of depth first searches
    for (int v : order)
        if (!mark[v]) {
            vector<int> component;
            scc_dfs(v, ag_rev, component);
            // printVector(component);
            int root = *min_element(begin(component), end(component));
            for (int u : component)
                aux_cc_comp[u] = root;
        }
    
    ag_rev.clear(); 
    // vector<vector<int>>().swap(ag_rev);
}

// ****************************************************************************************************************************************************** 

int main(int argc, char* argv[])
{   
    // ************************************
    // *** IO + data preparation ***
    // ************************************
    {
        ios_base::sync_with_stdio(false); cin.tie(0); cout.tie(0); // Fast IO

        CLI::App app{"panbubble"};
        auto decomp = app.add_subcommand("decompose", "Decompose GFA file into panbubbbles");
        // if(argc != 3){
        //     return 1;
        // }
        // inputpath = argv[1];
        // outputdir = argv[2];
        
        decomp->add_option("-i, --input", inputpath, "Input GFA")->required();
        decomp->add_option("-o, --output", outputdir, "Directory for saving the output files")->required();
        decomp->add_option("-d, --depth", maxdepth, "Maximum depth up to which the panbubbles are to be reported (1 means outermost)")->default_val(1);
        decomp->add_option("-m, --maxsize", maxsize, "Output vertices in the panbubble of size atmost [maxsize]")->default_val(1000);
        decomp->add_flag("-a, --allele", print_allele, "Whether alleles are to be reported for each panbubble");
        decomp->add_flag("-c, --cycle-equivalent", print_equivalent, "Whether cycle equivalent pairs are to be reported");
        decomp->add_flag("-r, --report-hairpins", print_hairpin, "Whether hairpins are to be reported");
        decomp->add_flag("-p, --print-panbubble-tree", print_panbubble_tree, "Whether the panbubble tree is to be printed");
         
        CLI11_PARSE(app, argc, argv);

        if(!fs::exists(outputdir))fs::create_directories(outputdir);
        
        summarypath = outputdir + "/summary.txt";
        freopen(summarypath.c_str(), "w", stdout);

        get_ne();

        cout << "Number of nodes in the input: " << n << endl;
        cout << "Number of edges in the input: " << edges << endl;

        // for a node x (0-indexed) in pangene graph - two nodes 2 * x (tail of arrow), 2 * x + 1 (head of arrow) are created in the bi-edged graph 
        // using vectors will be good coz we will have to add capping backedges as well
        g.resize(2 * n); // +delta is for S if required -> not needed (choose S as one of the tip ends)
        has_self_loop.resize(2 * n);
        // adding black edges first <- important since don't want black edge to appear as a back or front edge
        for(int i = 0; i < n; i++){
            int n1 = i << 1, n2 = n1 + 1;
            g[n1].pb({n2, -1}); g[n2].pb({n1, -1});
        }
        make_graph(); // adding grey edges
    }

    // ************************************
    // *** Initialisation ***    
    // ************************************
    {
        // Graph compression - removing unitigs --- Not implementing currently -> can get done without this
    
        // ************************************
        // *** Using Tarjan's algorithm to find bridge edges ***  // could have used sese but this is much faster in practice 
        // Implementation adapted from : https://cp-algorithms.com/graph/bridge-searching.html
        // ************************************
        {   
            // ** Initialise **
            mark.resize(2 * n);
            mark_bridge.resize(n);
            tin.resize(2 * n, -1); low.resize(2 * n, -1);

            for(int i = 0; i < 2 * n; i++){
                if(mark[i]){
                    continue;
                }
                dfs_bridge(i, -1);
            }
            
            // printArgs("Bridges Found:");
            // for(int i = 0; i < n; i++){
            //     if(mark_bridge[i])cout << i << " ";
            // }
            // cout << endl;

            // ** Clearing the memory **
            tin.clear(); low.clear();
            // vector<int>().swap(tin); vector<int>().swap(low);
        }

        // ************************************
        // *** Splitting the graph from the bridge point ***  // done carefully such that unitig compression is not required
        // ************************************
        {   
            // ************************************
            // *** Finding new n and marking the type of split ***
            // ************************************
            {
                // ** Initialise **
                n_processed = 2 * n;
                type_edge.resize(n);

                for(int i = 0; i < n; i++){ // iterate over all the black edges
                    if(mark_bridge[i]){
                        int u = i << 1, v = u + 1;
                        set<int> neigh_uv = find_unique_excluding_node(u, v); int cnt_uv = neigh_uv.size();
                        set<int> neigh_vu = find_unique_excluding_node(v, u); int cnt_vu = neigh_vu.size();
                        bool splitfrom_u = false, splitfrom_v = false;
                        if(cnt_uv >= 2){
                            splitfrom_u = true;
                        }else if(cnt_uv == 1){
                            int x = *neigh_uv.begin();
                            set<int> neigh_xu = find_unique_excluding_node(x, u); int cnt_xu = neigh_xu.size();
                            if(cnt_xu >= 2)splitfrom_u = true;
                        }
                        if(cnt_vu >= 2){
                            splitfrom_v = true;
                        }else if(cnt_vu == 1){
                            int w = *neigh_vu.begin();
                            set<int> neigh_wv = find_unique_excluding_node(w, v); int cnt_wv = neigh_wv.size();
                            if(cnt_wv >= 2)splitfrom_v = true;
                        }
                        if(splitfrom_u && splitfrom_v){
                            n_processed += 2;
                            type_edge[i] = '3';
                        }else if(splitfrom_u){
                            if(cnt_vu == 0){ // tip case
                                if(!has_self_loop[v]){
                                    type_edge[i] = '0';
                                }else{
                                    n_processed += 2;
                                    type_edge[i] = '3';
                                }
                            }else{
                                n_processed++;
                                type_edge[i] = '1';
                            }
                        }else if(splitfrom_v){
                            if(cnt_uv == 0){
                                if(!has_self_loop[u]){
                                    type_edge[i] = '0';
                                }else{
                                    n_processed += 2;
                                    type_edge[i] = '3';
                                }
                            }else{
                                n_processed++;
                                type_edge[i] = '2';
                            }                        
                        }else{
                            type_edge[i] = '4';
                        }
                    }else type_edge[i] = '0';

                    // cout << i << " " << type_edge[i] << endl;
                }

                cout << "Number of nodes after initialisation: " << (n_processed / 2) << endl;

                // printArgs("Type of edges in the original graph");
                // printVector(type_edge);
                
                // ** Clearing the memory **
                mark_bridge.clear(); has_self_loop.clear();
                // vector<bool>().swap(mark_bridge); vector<bool>().swap(has_self_loop);
            }

            // ************************************
            // *** Initialising the processed graph ***
            // ************************************
            {   
                // ** Initialise **
                id_ptr = 2 * n;
                // edges = 0;
                dual.resize(n_processed);
                g_processed.resize(n_processed);
                maxd = n_processed + 100;

                for(int i = 0; i < n; i++){ // copying data from graph g
                    int u = i << 1, v = u + 1;
                    switch (type_edge[i])
                    {
                        case '0':
                            neigh_copy(u, -1); neigh_copy(v, -1); break;
                        case '1':
                            neigh_copy(u, v); break;
                        case '2':
                            neigh_copy(v, u); break;
                        case '3':
                            neigh_copy(u, v); neigh_copy(v, u); break;
                        default:
                            break;
                    }
                }
                // assert(id_ptr == n_processed);
                // edges /= 2;
                // cout << "Number of edges after initialisation: " << edges << endl;

                // printGraph(g);
                // printGraph(g_processed);

                // ** Clearing the memory **
                g.clear(); type_edge.clear();
                // vector<vector<pii>>().swap(g); vector<char>().swap(type_edge);
            }
        }
        
        // ************************************
        // *** Finding number of connected components in the processed graph***    
        // count of vertices = n_processed
        // count of black edges = n (no of genes)
        // ************************************
        {
            // ** Initialise **
            mark.resize(n_processed);
            fill(mark.begin(), mark.end(), false);
            biedged_connected_comp.resize(n_processed);

            for(int i = 0; i < n_processed; i++){
                if(mark[i] || g_processed[i].size() == 0){
                    continue;
                }
                dfs_comp(i);
                first_node.pb(i);
                after_initialisation_comp_cnt++;
            }
            cout << "Number of components in the processed graph: " << after_initialisation_comp_cnt << endl;
        }

        // ************************************
        // *** Identifying tips ***
        // ************************************
        {
            // ** Initialise **
            tips.resize(after_initialisation_comp_cnt);
            
            for(int i = 0; i < n_processed; i++){
                if(g_processed[i].size() == 0)continue;
                // only connected via a black edge or a grey edge to (i ^ 1)
                if(find_unique_including_selfloop(i) == 0){ // much stricter condition than below
                // if(g_processed[i].size() == 1){
                    tips[biedged_connected_comp[i]].pb(i);
                    // assert((id_in_original_graph[g_processed[i][0].F] ^ id_in_original_graph[i]) == 1); // checking whether the node is connected via a black edge
                }
            }

            int cnt_zero = 0;
            for(int i = 0; i < after_initialisation_comp_cnt; i++){
                // printArgs("Tips", i);
                // printVector(tips[i]);
                if(tips[i].size() == 0)cnt_zero++;
            }
            cout << "Number of components with zero tips: " << cnt_zero << endl;
        }
    }

    // ************************************
    // *** Finding Possible Bibubble Pairs ***
    // ************************************
    {
        // ************************************
        // *** Finding node S ***
        // ************************************        
        {
            // ** Initialise **
            tip_start.resize(after_initialisation_comp_cnt, -1);
            Start.resize(after_initialisation_comp_cnt, -1);
            depth.resize(n_processed);
            bl.resize(n_processed);
            rm_cnt.resize(n_processed);
            fill(mark.begin(), mark.end(), false);
            
            for(int it = 0; it < after_initialisation_comp_cnt; it++){
                if(tips[it].size() == 0){ 
                    int val = edges;
                    sese_minbracket(first_node[it], -1, Start[it], val); // node belonging to the edge with minimum bracket set size
                    // printArgs("Cyclic component Starting point", ilmap[id_in_original_graph[Start[it]] >> 1]);
                    // no extra edges required
                }else{
                    Start[it] = tips[it][0]; 
                    // if(tips.size() > 1){// else no need to add an extra edge
                        for(int i = 1; i < tips[it].size(); i++){
                            if(tip_start[it] == -1)tip_start[it] = tot_grey; // marking the starting id for the edges added because of tips
                            add_edge(g_processed, Start[it], tips[it][i]);
                        }
                    // }
                }
                if(tips[it].size() <= 1)tip_start[it] = maxv; 
            }

            // printArgs("S nodes:");
            // printVector(Start);

            // printGraph(g_processed);

            // ** Clearing the memory **
            rm_cnt.clear(); first_node.clear();
            // vector<int>().swap(rm_cnt); vector<int>().swap(first_node);
        }
    
        // ************************************
        // *** Finding depth and number of backedges for individual components ***
        // ************************************
        {
            // ** Initialise **
            backedge_cnt.resize(after_initialisation_comp_cnt);
            fill(depth.begin(), depth.end(), 0);
            fill(mark.begin(), mark.end(), false);

            for(int it = 0; it < after_initialisation_comp_cnt; it++){   
                // Finding an upperbound on the number of backedges <- only these will contribute to the brackets
                nodes = 0;
                dfs_depth(Start[it], -1, it);
                backedge_cnt[it] += nodes; // capping backedges (for every node at most one capping back edge)
            }
        }

        // ************************************
        // *** Reinitialising edge ids ***
        // ************************************
        {   
            tot_grey_comp.resize(after_initialisation_comp_cnt);

            for(int pass = 0; pass < 2; pass++){
                fill(mark.begin(), mark.end(), false);

                for(int it = 0; it < after_initialisation_comp_cnt; it++){
                    if(pass == 1 && tip_start[it] == maxv)continue; // no change required
                    reinitialise_edge_ids(Start[it], -1, it, pass);
                }
            }

            tot_grey_comp.clear();
        }

        // ************************************
        // *** SESE ***
        // ************************************
        {
            remove_brackets.resize(n_processed);
            // visit_time.resize(n_processed, -1);
            bl.clear();
            fill(mark.begin(), mark.end(), false);
            for(int it = 0; it < after_initialisation_comp_cnt; it++){
                multiplier = 1;
                while(multiplier <= backedge_cnt[it]){
                    multiplier *= 10;
                }
                st.clear(); 
                // st0 = -1;
                sese(Start[it], -1); // graph g is not changed
                // assert(stack_trace.size() == 0); 

                // corner case
                if(tips[it].size() == 1){
                    canonical_sese.pb({dual[tips[it][0]], dual[tips[it][0]]});
                }
            }
    
            // ** Clearing the memory **
            bl.clear(); remove_brackets.clear(); // visit_time.clear(); 
            Start.clear(); depth.clear(); 
            backedge_cnt.clear(); tips.clear();
            // vector<bracketlist*>().swap(bl); vector<vector<edge*>>().swap(remove_brackets); vector<int>().swap(Start); vector<int>().swap(depth);  
            // vector<int>().swap(backedge_cnt);  vector<vector<int>>().swap(tips);
        }
        
        // ************************************
        // *** Printing result ***
        // ************************************
        {
            possible_pairs = canonical_sese.size();
            if(print_equivalent){
                summarypath = outputdir + "/cycle_equivalent.txt";
                freopen(summarypath.c_str(), "w", stdout);

                cout << "Total cycle equivalent pairs found: " << possible_pairs << endl;
                
                for(int i = 0; i < possible_pairs; i++){
                    pii rs = canonical_sese[i];
                    cout << get_label(rs.F, rs.S) << endl;                  
                }
            }
        }
    }

    // ************************************
    // *** Finding Valid Bibubble Pairs in G_O ***
    // ************************************
    {   
        if(possible_pairs != 0){
            // ************************************
            // *** SCC ***
            // ************************************
            {   
                aux_cc_comp.resize(n_processed, -1);
                ag.resize(n_processed);
                // need not connect tips here
                make_auxillary_graph_from_biedged_graph(); // need not read the file again
                // printArgs("Auxiliary graph:");
                // printGraph(ag);
                scc();
            }

            // ************************************
            // *** Finding U set  (Using biedged graph representation for the same) ***
            // ************************************
            {
                fill(mark.begin(), mark.end(), false);
                bb_comp.resize(n_processed, -1);
                bb_nodes.resize(possible_pairs);
                considered.resize(n_processed, false);
                opp_entrance.resize(n_processed, -1);
                
                for(int i = 0; i < possible_pairs; i++){ // space requirement is linear since the inner-most nested bubbles will appear on the top
                    pii rs = canonical_sese[i];
                    // upd_node(id_in_original_graph[rs.F] ^ 1, i); // will create issues when two panbubbles are adjacent
                    mark_bb_nodes(rs.F, rs, i); // do not add the end points, checking only in one direction
                    mark[rs.F] = mark[rs.S] = mark[dual[rs.F]] = mark[dual[rs.S]] = true;
                    bb_comp[rs.F] = bb_comp[rs.S] = bb_comp[dual[rs.F]] = bb_comp[dual[rs.S]] = i;
                    opp_entrance[rs.F] = rs.S; opp_entrance[rs.S] = rs.F;
                }

                // for(int i = 0; i < possible_pairs; i++){
                //     printArgs("bb_nodes", i, ":");
                //     printVector(bb_nodes[i]);
                // }

                considered.clear(); opp_entrance.clear();
            }
        }

        valid.resize(possible_pairs, true);
        bg.resize(possible_pairs);
        root_bg.resize(possible_pairs, true);
        inb.resize(possible_pairs, false);
        valid_pairs = possible_pairs;

        // ************************************
        // *** Validating the Theorem ***
        // ************************************
        for(int i = 0; i < possible_pairs; i++){
            pii rs = canonical_sese[i];
            for(int u : bb_nodes[i]){
                // assert(bb_comp[u] != -1);
                // assert(!end_gene(u, rs));
                // if(!end_gene(u, rs)){ // check only for internal nodes
                    if(valid[bb_comp[u]] && !inb[bb_comp[u]]){
                        // cout << i << " " << u << " " << get_label(u) << endl;
                        if(end_gene(u, canonical_sese[bb_comp[u]])){
                            root_bg[bb_comp[u]] = false;
                            bg[i].pb(bb_comp[u]);
                            inb[bb_comp[u]] = true;
                        }
                    }
                    if(!((aux_cc_comp[u] == aux_cc_comp[rs.F] || aux_cc_comp[u] == aux_cc_comp[dual[rs.F]]) && valid[bb_comp[u]])){
                        valid[i] = false;
                        valid_pairs--;
                        break;
                    }
                // }
            }

            // Can contain hairpins with size = 0 (edge with a loop in the end)
            if(bb_nodes[i].size() == 0){
                // valid[i] = false;
                // valid_pairs--;
            }

            if(rs.F == rs.S)hairpins += valid[i];
        }

        bb_nodes.clear(); bb_comp.clear(); inb.clear();

        // ************************************
        // *** Printing the result ***
        // ************************************        
        {
            // ** Clearing the memory **
            tip_start.clear(); 
            // vector<ll>().swap(tip_start);

            if(print_panbubble_tree){
                summarypath = outputdir + "/panbubble_tree.txt";
                freopen(summarypath.c_str(), "w", stdout);
                for(int i = 0; i < possible_pairs; i++){
                    pii rs = canonical_sese[i];
                    // if(valid[i] && rs.F != rs.S){
                    if(valid[i]){
                        if(root_bg[i]){
                            cout << "_, " << get_label_from_id(i) << endl;
                        }
                        for(int v : bg[i]){
                            cout << get_label_from_id(i) << ", " << get_label_from_id(v) << endl;
                        }
                    }
                }
            }

            // ** Removing the bubbles with depth > maxdepth **
            bubble_depth.resize(possible_pairs, -1);
            for(int i = 0; i < possible_pairs; i++){
                pii rs = canonical_sese[i];
                if(valid[i] && rs.F != rs.S && root_bg[i]){
                    qbg.push_back(i);
                    bubble_depth[i] = 1;
                }
            }

            while(!qbg.empty()){
                int u = qbg.front();
                qbg.pop_front();
                for(int v : bg[u]){
                    if(bubble_depth[v] != -1 || !valid[v]){
                        continue;
                    }
                    bubble_depth[v] = bubble_depth[u] + 1;
                    qbg.push_back(v);
                }
            }

            for(int i = 0; i < possible_pairs; i++){
                if(bubble_depth[i] > maxdepth){
                    valid[i] = false;
                    valid_pairs--;
                }
            }

            summarypath = outputdir + "/summary.txt";
            freopen(summarypath.c_str(), "a", stdout);
            cout << "Total panbubbles found till the input depth: " << (valid_pairs - hairpins) << endl;

            summarypath = outputdir + "/panbubble.txt";
            freopen(summarypath.c_str(), "w", stdout);

            // string header = "CC	BB  bbID  parID  side1  side2  #alleles  #genes  geneList\nCC   AL  #hap  walk haplotypeList\nCC";
            // cout << header << endl;
            
            for(int i = 0; i < possible_pairs; i++){
                pii rs = canonical_sese[i];
                if(!valid[i] || rs.F == rs.S) continue;
                cout << get_label(rs.F, rs.S) << endl;
            }

            if(print_hairpin){
                summarypath = outputdir + "/summary.txt";
                freopen(summarypath.c_str(), "a", stdout);
                cout << "Total hairpins found: " << hairpins << endl;
    
                summarypath = outputdir + "/hairpins.txt";
                freopen(summarypath.c_str(), "w", stdout);

                for(int i = 0; i < possible_pairs; i++){
                    pii rs = canonical_sese[i];
                    if(!valid[i] || rs.F != rs.S) continue;
                    cout << get_label(rs.F, rs.S) << endl;
                }
            }
        }
    }
} 

// TODO:
// cyclic cases