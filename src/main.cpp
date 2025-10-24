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
const int maxv = INT32_MAX;
int tot_gray = 0; // for storing total number of gray edges, will help in accessing the brackets
int maxd; // max depth possible = total number of nodes
string summarypath; // path where the given stats will be written
vector<bool> mark; // for marking the vertices that are have been visited
vector<string> ilmap; // for storing the gene for a particular label
// map<int, int> id_in_original_graph; // used to find the id in original graph, in order to use the property of black edges map[u] ^ map[v] = 1

string get_single_label(int x, int ty){
    return ((ty == 0 ? (x & 1 ? ">" : "<") : (x & 1 ? "<" : ">")) + ilmap[x >> 1]);
}

string get_label(int x, int y){
    return "( " + (get_single_label(x, 0) + " " + get_single_label(y, 1)) + " )"; 
}
// ****************************************************************************************************************************************************** 

// ****************************************************************************************************************************************************** 
// **** Structures ****
// ****************************************************************************************************************************************************** 
struct edge{
    int id; // a bracket can uniquely be identified by the lower and higher vertices it connects to : low -> lower height, assigning a unique to each
    // int node; // required node corresponding to the black edge when edge 'id' is the top in the bracketlist
    // int sz; // required size of bracketlist corresponding to the black edge when edge 'id' is the top in the bracketlist
    edge* front;
    edge* back; // for doubly linked list implementation -> helps in O(1) deletion
    // edge(int eid):id(eid), node(-1), sz(-1), front(nullptr), back(nullptr){}
    edge(int eid):id(eid), front(nullptr), back(nullptr){}
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
// **** Helper functions ****
// ****************************************************************************************************************************************************** 
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

void printVector(vector<int>& v, int ty = 0){
    if(!PRINT)return;
    for(int x : v){
        cout << (ty == 0 ? to_string(x) : get_single_label(x, 0)) << " ";
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
// **** Input ****
// ****************************************************************************************************************************************************** 
fstream f;
string inputpath, outputdir;
bool print_allele = false, print_equivalent = false, print_hairpin = false, print_panbubble_tree = false; // whether to report alleles, cycle equivalent pairs and hairpins
int maxsize = 1000, maxdepth = 1; // max size of panbubbles whose vertices are also to be printed, max depth to which the panbubbles are to be reported
int n, edges; // no of nodes (genes), no of edges
// int before_initialisation_comp_cnt = 0; // no of connected components in the input graph before initialisation
vector<bool> has_self_loop; // whether the node has a self loop <- can't be cycle equivalent with any other node
vector<vector<pii>> g; // (node_id, gray_edge_id)
map<string, int> lmap; // for storing the label for a particular gene

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
    g[id1].pb({id2, tot_gray}); 
    if(id1 != id2) g[id2].pb({id1, tot_gray}); // 0-indexed
    tot_gray++;
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
// **** Compacted Graph ****
// ****************************************************************************************************************************************************** 
int cnt_gray_edge = 0; // count of gray edges in the compacted graph
int component = 0; // no of connected components in the compacted graph
vector<bool> can_compact; // whether a gene has been visited during clustering used for compaction
vector<bool> upd_vertex; // whether the adjacency list of a vertex has been updated during compaction
vector<int> rm_cnt; // vector storing the number of brackets ending at that node
vector<int> biedged_connected_comp; // connected-component in the compacted graph to which a node belongs to 
vector<int> tip_start; // will store the id from which extra edges because of tips are added for each component in the compacted graph
vector<int> first_node; // for storing a node in every connected-component in the compacted graph
vector<int> Start; // starting point of search for each component in compacted graph
vector<int> depth; // for storing the depth of the vertices in the spanning tree
vector<vector<int>> tips; // vector for storing the vertices that represent the tip
vector<vector<pii>> g_compacted; // (node_id, gray_edge_id)
queue<int> q_compact; // queue used for clustering the edges in a linear chain

bool chk_for_compaction(int id){
    int n1 = id << 1, n2 = n1 + 1;
    if(g[n1].size() == 2 && g[n2].size() == 2 && !has_self_loop[n1] && !has_self_loop[n2]){
        // cout << n1 << " " << n2 << endl;
        if((g[n1][0].F != g[n1][1].F) && (g[n2][0].F != g[n2][1].F)){
            int u = g[n1][1].F, v = g[n2][1].F;
            if(g[u].size() > 2 || g[v].size() > 2)return false;
            can_compact[id] = true;
            return true;
        }else return false;
    }else return false;
}

void dfs_comp(int u){ // Finding connected components in the compacted graph
    mark[u] = true;
    biedged_connected_comp[u] = component;
    for(pii child : g_compacted[u]){
        int v = child.F;
        if(mark[v])continue;
        dfs_comp(v);
    }
}

int find_unique_including_selfloop(int u){ // skip u ^ 1
    set<int> s;
    for(pii child : g_compacted[u]){
        if((child.F ^ u) == 1)continue;
        s.insert(child.F);
        break;
    }
    return s.size();
}

set<int> find_unique_excluding_node(int u, int v){// exclude node u and v
    set<int> s;
    for(pii child : g_compacted[u]){
        if(child.F == u || child.F == v)continue;
        s.insert(child.F);
        if(s.size() == 2)break;
    }
    return s;
}
// ****************************************************************************************************************************************************** 

// ****************************************************************************************************************************************************** 
// **** Single Entry Single Exit Algorithm (Johnson et al., 1994) ****
// ******************************************************************************************************************************************************  
// int nodes = 0; // no of nodes in a given component in the input graph
// int st0 = 0; // top node when the bracketlist is empty
int time = 0;
ll multiplier = 1; // for computing the hash
int number_of_classes = 0; // storing the number of valid panbubble pairs

vector<int> entry_time, exit_time; // for storing the entry and exit time of vertices during dfs
vector<int> stack_trace; // for storing the vertices in the stack during dfs traversal
vector<pii> possible_hairpins; // for storing possible hairpin's entrance vertex
// vector<int> bubble_depth; // depth of the bubbles in the bubble tree
// vector<int> bracket_seq; // for storing the id of the bracket sequences so as to make the bubble tree
// vector<int> visit_time; // for storing the visit time of the vertices during the dfs traversal
// vector<bool> root_bg; // root bubble id's in the bubble tree
// vector<bool> inb; // whether the node has been added in the bubble tree
vector<bool> is_bridge; // whether the edge is a bridge in G_U
map<pii, vector<pii>> canonical_sese; // for storing the canonical sese pairs in a given class
// vector<vector<int>> bg; // vector for storing the adjacency matrix for the bubble graph
// map<ll, int> st; // (edge, size) -> hashed into a ll key
map<pii, int> st; // (id, size) -> as the key
vector<bracketlist*> bl; // vector storing the bracket list for different nodes
vector<vector<edge*>> remove_brackets; // vector that stores the brackets to be removed while exiting a particular vertex during SESE algorithm
// deque<int> qbg; // queue maintained while doing bfs on the bubble tree

// ll get_key(int id, int size){
//     return size * multiplier + id;
// }

int find_unique_excluding_selfloop(int u, int v){// exclude nodes u, u ^ 1, v
    set<int> s;
    for(pii child : g_compacted[u]){
        if((child.F ^ u) <= 1 || child.F == v)continue;
        s.insert(child.F);
        break;
    }
    return s.size();
}

void dfs_depth(int u, int parent){
    // assert(!mark[u]);
    if(parent != -1)depth[u] = depth[parent] + 1;
    mark[u] = true;
    for(pii child : g_compacted[u]){
        int v = child.F;
        if(mark[v]){
            // if(depth[v] < depth[u] && v != parent)backedge_cnt[comp]++;
            continue;
        }
        dfs_depth(v, u);
    }
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

void sese_minbracket(int u, int parent, int& x, int& val){
    // printArgs(u, parent, x, val);
    // assert(!mark[u]);
    mark[u] = true;
    if(parent != -1)depth[u] = depth[parent] + 1;

    for(pii child : g_compacted[u]){
        int v = child.F;
        if(v == parent || v == u)continue; // v = u won't contribute anything and v == parent won't contribute to backedge
        if(mark[v]){// back-edge -> will come first in the dfs traversal for the node at greater depth
            if(depth[v] > depth[u])continue; // front-edge -- multi edges will be taken care of here
            // assert(child.S != -1); // black edges can't be backedges
        }else sese_minbracket(v, u, x, val);
    }

    bl[u] = new bracketlist(); 

    for(pii child : g_compacted[u]){
        int v = child.F;
        if(depth[v] != depth[u] + 1)continue;
        if(bl[v] && !bl[v] -> merged){
            // printArgs(u, v, bl[v]->sz);
            if((u ^ v) == 1){
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
    for(pii child : g_compacted[u]){
        int v = child.F;
        if(depth[v] >= depth[u] - 1)continue;
        rm_cnt[v]++;
        bl[u]->sz++;
    }
}

void sese(int u, int parent){
    stack_trace.pb(u);
    mark[u] = true;
    entry_time[u] = time++;

    bool can_end = true;
    int h0 = maxd, h1 = maxd, h2 = maxd;

    for(pii child : g_compacted[u]){
        int v = child.F;
        if((v ^ u) == 1 && child.S != -1){ // a parallel gray edge <- can not be the panbubble end
            can_end = false;
            break;
        }
    }

    for(pii child : g_compacted[u]){
        int v = child.F;
        if(v == parent || v == u)continue; // v = u won't contribute anything and v == parent won't contribute to backedge
        if(mark[v]){// back-edge -> will come first in the dfs traversal for the node at greater depth
            if(depth[v] > depth[u])continue; // front-edge -- multi edges will be taken care of here
            // assert(child.S != -1); // black edges can't be backedges
            h0 = min(h0, depth[v]);
        }else sese(v, u);
    }

    bl[u] = new bracketlist(); 

    for(pii child : g_compacted[u]){
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
    for(pii child : g_compacted[u]){
        int v = child.F;
        if(depth[v] >= depth[u] - 1)continue;
        edge* ed = new edge(child.S); // tot_gray count starts from 0
        remove_brackets[v].pb(ed);
        merge(bl[u], new bracketlist(1, depth[v], ed, ed));
    }

    // capping back edge
    if(h2 < h0){
        // printArgs("finding capping backedge");
        int w = stack_trace[h2];
        edge* ed = new edge(cnt_gray_edge); cnt_gray_edge++; // need not modify g_compacted, will not be traversing along this edge
        remove_brackets[w].pb(ed);
        merge(bl[u], new bracketlist(1, depth[w], ed, ed));
        // printArgs("capping back edge:", get_single_label(u, 0), get_single_label(w, 0));
    }

    // finding cycle equivalence
    if(parent != -1 && ((parent ^ u) == 1) && can_end){
        // printArgs("finding cycle equivalence:", u, parent);
        // ll key;
        pii key;
        if(bl[u]->sz > 0){
            int br = bl[u]->end->id;
            // assert(br != -1); // not a black edge
            // key = get_key(br, bl[u]->sz);
            key = {br, bl[u]->sz};
        }else{
            // key = 0;
            key = {0, 0};
        }

        if(st.find(key) != st.end()){
            int w = st[key];
            if(find_unique_excluding_selfloop(u, w) == 1 && find_unique_excluding_selfloop(w, u) == 1){
                // printArgs("cycle equivalent pair:", u, w);
                canonical_sese[key].pb({u, w});
            }
        }
        st[key] = parent; // will happen regardless you found something or not
    }

    if(bl[u]-> sz == 0){
        is_bridge[u >> 1] = true;
    }

    exit_time[u] = time++;
    stack_trace.rb();
}
// ******************************************************************************************************************************************************  

// ****************************************************************************************************************************************************** 
// **** Brute force ****
// ******************************************************************************************************************************************************  
int counter = 0; // counter used in the brute force checking
vector<int> considered; // for saving the vertices visited during the brute force check for a panbubble
vector<short> to_chk; // which side of gene to check during second pass of brute (-1 -> initialisation, 0 -> 2 * i, 1 -> 2 * i + 1, 2 -> both)
vector<short> visit_order; // from which side the gene has been visited (-1 -> not visited, 0 -> rs.F, 1 -> rs.S, 2 -> both)

bool end_gene(int& u, pii& rs){
    return u == rs.F || u == (rs.F ^ 1) || u == rs.S || u == (rs.S ^ 1);
}

bool root_leaf_path(int u, int v){
    return entry_time[u] <= entry_time[v] && exit_time[u] >= exit_time[v];
}

void upd_node(int u, int id, int ty){
    visit_order[u] = ty;

    int gene = u >> 1;    
    if(ty == 0){
        if(to_chk[gene] == -1){
            counter++;
            to_chk[gene] = u & 1;
            considered.pb(u);
        }else{
            if((u & 1) != to_chk[gene]){
                to_chk[gene] = 2;
            }
        }
    }else{
        if(to_chk[gene] == 2){
            counter--;
            to_chk[gene] = -2; // -2 -> same gene does not reduce the counter multiple times
        }else{
            short required = 1 - (u & 1);
            if(to_chk[gene] == required || to_chk[gene] == -1){
                counter--;
                to_chk[gene] = -2;
                if(to_chk[gene] == -1)
                    considered.pb(u);
            }
        }
    }
}

bool mark_nodes(int u, pii& rs, int id, int ty){
    if(is_bridge[u >> 1] && (!root_leaf_path(rs.F, u) || !root_leaf_path(u, rs.S)))return false;

    if(!end_gene(u, rs)){
        upd_node(u, id, ty);
    }

    for(pii child : g_compacted[u ^ 1]){
        int v = child.F;

        if(child.S >= tip_start[biedged_connected_comp[rs.S]] || child.S == -1 || end_gene(v, rs))continue; // don't consider the edges added because of tips
        if(visit_order[v] == ty || (ty == 1 && visit_order[v] == 0)){ // second case to take care of the hairpins
            continue;
        }
        if(!mark_nodes(v, rs, id, ty))return false;
    }
    return true;
}
// ****************************************************************************************************************************************************** 

// ****************************************************************************************************************************************************** 
// **** Output ****
// ****************************************************************************************************************************************************** 
vector<pii> valid_panbubbles, valid_hairpins; // for storing valid panbubbles and hairpins
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
        decomp->add_flag("-c, --cycle-equivalent", print_equivalent, "Whether cycle equivalent classes are to be reported");
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
        maxd = 2 * n + 100;

        // adding black edges first <- important since don't want black edge to appear as a back or front edge
        for(int i = 0; i < n; i++){
            int n1 = i << 1, n2 = n1 + 1;
            g[n1].pb({n2, -1}); g[n2].pb({n1, -1});
        }
        make_graph(); // adding gray edges
    }

    // ************************************
    // *** Graph Cleaning ***    
    // ************************************
    {
        // ************************************
        // *** Removing parallel gray edges ***
        // ************************************
        {
            for(int i = 0; i < 2 * n; i++){
                sort(g[i].begin() + 1, g[i].end());

                // case for a black and a gray edge
                int add = 0;
                if(g[i].size() >= 2){
                    if(g[i][0].F == g[i][1].F){
                        add = 1;
                    }
                }
                g[i].erase(unique(g[i].begin() + add, g[i].end(), [](auto &a, auto &b){return a.first == b.first;}), g[i].end());
            }
        }

        // ************************************
        // *** Graph compaction - removing unitigs (removing based on the genes (edges) and not the vertices) ***
        // ************************************
        {   
            // clustering till the end edges are achieved
            can_compact.resize(n);
            upd_vertex.resize(2 * n);
            g_compacted.resize(2 * n);

            for(int i = 0; i < n; i++){
                if(can_compact[i])continue;
   
                if(chk_for_compaction(i)){
                    q_compact.push(i);
                    int left = i << 1, right = left + 1;
                    while(!q_compact.empty()){
                        int gene = q_compact.front(); q_compact.pop();

                        int gene_L = gene << 1, gene_R = gene_L + 1;
                        int add_L = g[gene_L][1].F >> 1, add_R = g[gene_R][1].F >> 1;
                        if(!can_compact[add_L]){
                            left = g[gene_L][1].F;
                            if(chk_for_compaction(add_L)){
                                q_compact.push(add_L);
                            }                            
                        }
                        if(!can_compact[add_R]){
                            right = g[gene_R][1].F;
                            if(chk_for_compaction(add_R)){
                                q_compact.push(add_R);
                            }
                        }
                    }

                    int l1 = left, l2 = left ^ 1;
                    g_compacted[l1].pb({l2, -1}); g_compacted[l2].pb({l1, -1});

                    int r1 = right, r2 = right ^ 1;
                    if((l1 ^ r1) != 1){
                        g_compacted[r1].pb({r2, -1}); g_compacted[r2].pb({r1, -1});

                        g_compacted[l1].pb({r1, cnt_gray_edge}); g_compacted[r1].pb({l1, cnt_gray_edge++});
                    }
                    upd_vertex[l1] = upd_vertex[r1] = true;

                    for(int j = 1; j < g[l2].size(); j++){ // have removed parallel gray edges
                        int u = g[l2][j].F;
                        if(upd_vertex[u] || (can_compact[u >> 1] && (u != r2)))continue; // u != r2 -> gray edge (cyclic case)
                        g_compacted[l2].pb({u, cnt_gray_edge});
                        g_compacted[u].pb({l2, cnt_gray_edge++}); // 0th edge maynot be the -1 edge
                    }
                    for(int j = 1; j < g[r2].size(); j++){
                        int u = g[r2][j].F;
                        if(upd_vertex[u] || can_compact[u >> 1])continue;
                        g_compacted[r2].pb({u, cnt_gray_edge});
                        g_compacted[u].pb({r2, cnt_gray_edge++});
                    }
                    upd_vertex[l2] = upd_vertex[r2] = true;

                    can_compact[l1 >> 1] = can_compact[r1 >> 1] = true;
                }
            }

            for(int i = 0; i < n; i++){
                if(can_compact[i])continue;

                int left = i << 1, right = left + 1;
                g_compacted[left].pb({right, -1}); g_compacted[right].pb({left, -1});
                for(int j = 1; j < g[left].size(); j++){
                    int u = g[left][j].F;
                    if(upd_vertex[u] || can_compact[u >> 1])continue;
                    g_compacted[left].pb({u, cnt_gray_edge}); 
                    g_compacted[u].pb({left, cnt_gray_edge++}); 
                }
                upd_vertex[left] = true;

                for(int j = 1; j < g[right].size(); j++){
                    int u = g[right][j].F;
                    if(upd_vertex[u] || can_compact[u >> 1])continue;
                    g_compacted[right].pb({u, cnt_gray_edge}); 
                    g_compacted[u].pb({right, cnt_gray_edge++});
                }
                upd_vertex[right] = true;
            }

            cout << "Number of gray edges after compaction: " << cnt_gray_edge << endl;

            // swapping the black edge so that it is the 0-indexed edge
            for(int i = 0; i < 2 * n; i++){
                for(pii& child : g_compacted[i]){
                    if(child.S == -1){
                        swap(child, g_compacted[i][0]);
                        break;
                    }
                }
            }

            // printGraph(g);
            // printArgs("------");
            // printGraph(g_compacted);

            // ** Clearing the memory **
            g.clear();
        }
    }
        
    // ************************************
    // *** Finding number of connected components in the compacted graph ***    
    // count of vertices = 2 * n (remains the same, although some vertices will have |adj[vertex]| = 0)
    // ************************************
    {
        // ** Initialise **
        mark.resize(2 * n);
        fill(mark.begin(), mark.end(), false);
        biedged_connected_comp.resize(2 * n);

        for(int i = 0; i < 2 * n; i++){
            if(mark[i] || g_compacted[i].size() == 0){
                continue;
            }
            dfs_comp(i);
            first_node.pb(i);
            component++;
        }
        cout << "Number of components in the compacted graph: " << component << endl;
    }

    // ************************************
    // *** Identifying tips ***
    // ************************************
    {
        // ** Initialise **
        tips.resize(component);
        
        for(int i = 0; i < 2 * n; i++){
            if(g_compacted[i].size() == 0)continue;
            // only connected via a black edge or a gray edge to (i ^ 1)
            if(find_unique_including_selfloop(i) == 0){ // much stricter condition than below
            // if(g_compacted[i].size() == 1){
                tips[biedged_connected_comp[i]].pb(i);
            }
        }

        int cnt_zero = 0;
        for(int i = 0; i < component; i++){
            // printArgs("Tips", i);
            // printVector(tips[i], 1);
            if(tips[i].size() == 0)cnt_zero++;
        }
        cout << "Number of components with zero tips: " << cnt_zero << endl;
    }

    // ************************************
    // *** Finding edges belonging to individual cycle equivalent classes ***
    // ************************************
    {
        // ************************************
        // *** Finding node S and adding extra gray edges ***
        // ************************************        
        {
            // ** Initialise **
            tip_start.resize(component, -1);
            Start.resize(component, -1);
            depth.resize(2 * n);
            bl.resize(2 * n);
            rm_cnt.resize(2 * n);
            fill(mark.begin(), mark.end(), false);
            
            for(int it = 0; it < component; it++){
                if(tips[it].size() == 0){ // recheck for EC22
                    int val = edges;
                    sese_minbracket(first_node[it], -1, Start[it], val); // node belonging to the edge with minimum bracket set size
                    // printArgs("Cyclic component Starting point", ilmap[Start[it] >> 1]);
                    // no extra edges required
                }else{
                    Start[it] = tips[it][0]; 
                    // if(tips.size() > 1){// else no need to add an extra edge
                        int u = Start[it];
                        for(int i = 1; i < tips[it].size(); i++){
                            if(tip_start[it] == -1)tip_start[it] = cnt_gray_edge; // marking the starting id for the edges added because of tips
                            int v = tips[it][i];
                            g_compacted[u].pb({v, cnt_gray_edge}); g_compacted[v].pb({u, cnt_gray_edge++});
                        }
                    // }
                }
                if(tips[it].size() <= 1)tip_start[it] = maxv; 
            }

            // printArgs("S nodes:");
            // printVector(Start);

            // ** Clearing the memory **
            rm_cnt.clear(); first_node.clear();
        }
    
        // ************************************
        // *** Finding depth for individual components ***
        // ************************************
        {
            // ** Initialise **
            fill(depth.begin(), depth.end(), 0);
            fill(mark.begin(), mark.end(), false);
            entry_time.resize(2 * n); exit_time.resize(2 * n);

            for(int it = 0; it < component; it++){   
                // Finding an upperbound on the number of backedges <- only these will contribute to the brackets
                dfs_depth(Start[it], -1);
            }
        }

        // ************************************
        // *** SESE ***
        // ************************************
        {
            remove_brackets.resize(2 * n);
            is_bridge.resize(n);
            // visit_time.resize(n_processed, -1);
            bl.clear();
            fill(mark.begin(), mark.end(), false);
            // int backedge_cnt = cnt_gray_edge + 2 * n; // upper bound
            
            // multiplier = 1;
            // while(multiplier <= backedge_cnt){
            //     multiplier *= 10;
            // }
            
            for(int it = 0; it < component; it++){
                st.clear(); 
                // st0 = -1;
                sese(Start[it], -1); // graph g is not changed
                
                // corner case
                if(tips[it].size() == 1){
                    possible_hairpins.pb({tips[it][0]^1, tips[it][0]^1});
                }
            }
            
            // printVector(entry_time); printVector(exit_time); printVector(is_bridge);

            // ** Clearing the memory **
            bl.clear(); remove_brackets.clear(); // visit_time.clear(); 
            Start.clear(); depth.clear(); 
            tips.clear();
        }
        
        // ************************************
        // *** Printing result ***
        // ************************************
        {
            number_of_classes = canonical_sese.size();
            if(print_equivalent){
                summarypath = outputdir + "/cycle_equivalent_classes.txt";
                freopen(summarypath.c_str(), "w", stdout);

                cout << "Total cycle equivalent classes found: " << number_of_classes << endl;
                
                for(const auto &[key, vec] : canonical_sese){
                    for(const pii &rs : vec){
                        cout << get_label(rs.F, rs.S) << " "; 
                    }
                    cout << endl;
                }
            }
        }
    }

    // ************************************
    // *** Finding valid panbubbles ***
    // ************************************
    {   
        fill(mark.begin(), mark.end(), false);
        visit_order.resize(2 * n, -1);
        to_chk.resize(n, -1);

        for(const auto &[key, vec] : canonical_sese){
            int sz = vec.size();
            for(int i = 0; i < sz; ){
                bool is_valid = true;
                int j;
                considered.clear();

                for(j = i; j < sz; j++){
                    pii rs = {vec[j].F, vec[i].S};

                    if(!root_leaf_path(rs.F, rs.S))continue;
                    
                    counter = 0;
                    
                    is_valid = mark_nodes(rs.F ^ 1, rs, i, 0); // do not add the end points, checking only in one direction
                    is_valid &= mark_nodes(rs.S ^ 1, rs, i, 1);
                    is_valid &= (counter == 0);

                    if(is_valid){
                        valid_panbubbles.pb(rs);
                        i = j + (j == i ? 1 : 0); 
                        break;
                    }
                }
                if(!is_valid || j == sz)i++;

                // ** Reinitialisation **
                for(int u : considered){
                    visit_order[u] = visit_order[u ^ 1] = -1;
                    to_chk[u >> 1] = -1;
                }
            }
        }

        // ** Clearing the memory **
        considered.clear(); visit_order.clear(); to_chk.clear(); 
    }
} 

// TODO:
// cyclic cases
// sending vars by refs
// clearing the memory