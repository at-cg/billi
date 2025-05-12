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
    edge* front;
    edge* back; // for doubly linked list implementation -> helps in O(1) deletion
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

void printVector(vector<int>& v){
    if(!PRINT)return;
    for(int x : v){
        cout << x << " ";
    }
    cout << endl;
}

void printVector(vector<pii>& v){
    if(!PRINT)return;
    for(pii x : v){
        cout << "(" << x.F << ", " << x.S << ") ";
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
        cout << i << ": ";
        printVector(g[i]);
    }
}
// ****************************************************************************************************************************************************** 

// ****************************************************************************************************************************************************** 
// **** Input ****
// ****************************************************************************************************************************************************** 
fstream f;
string inputpath, outputdir;
int n, edges; // no of nodes (genes), no of edges
// int before_initialisation_comp_cnt = 0; // no of connected components in the input graph before initialisation
vector<string> ilmap; // for storing the gene for a particular label
vector<bool> has_self_loop; // whether the node has a self loop <- can't be cycle equivalent with any other node
// vector<pii> eid; // stores the edges
vector<vector<pii>> g; // (node_id, grey_edge_id)
map<string, int> lmap; // for storing the label for a particular gene

void get_ne(){
    f.open(inputpath, ios::in);
    string line;
    regex strip("^\\s+|\\s+$"), split("\\t");
    vector<string> tokens;
    n = 0; edges = 0;
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
                // assert(tokens.size() >= 3);
                lmap[tokens[1]] = n; // 0-indexed
                ilmap.pb(tokens[1]);
                n++;
            }else if(tokens[0] == "L"){
                // assert(tokens.size() >= 5);
                edges++;
            }
        }
        f.close();
    }
}

void add_edge(vector<vector<pii>>& g, int& id1, int& id2){
    g[id1].pb({id2, tot_grey}); 
    if(id1 != id2) g[id2].pb({id1, tot_grey}); // 0-indexed
    // eid.pb({min(id1, id2), max(id1, id2)});
    tot_grey++;
}

// void printEdges(){
//     for(int i = 0; i < eid.size(); i++){
//         printArgs(i, eid[i].F, eid[i].S);
//     }
// }

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
                // assert(tokens.size() >= 5);
                int n1 = lmap[tokens[1]]; string s1 = tokens[2]; 
                int n2 = lmap[tokens[3]]; string s2 = tokens[4]; 
                // if(n1 == n2 && s1 == s2)continue; // won't use this edge for traversal

                // making sure edges are added only in one direction <- consistency of sign
                if(n1 > n2){
                    // continue;
                    swap(n1, n2); 
                    if(s1 == "+")s1 = "-";
                    else s1 = "+";
                    if(s2 == "+")s2 = "-";
                    else s2 = "+";
                    swap(s1, s2);
                }

                // n1 and n2 are 0-indexed
                int id1 = n1 << 1, id2 = n2 << 1;
                if(s1 == "+")id1++;
                if(s2 == "-")id2++;
                add_edge(g, id1, id2);
                // printArgs(id1, id2, tot_grey);

                if(n1 == n2 && s1 != s2){
                    // assert(id1 == id2);
                    has_self_loop[id1] = true;
                }
            }
        }
        f.close();
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
vector<vector<pii>> g_processed; // (node_id, grey_edge_id)
map<int, int> id_in_original_graph; // used to find the id in original graph, in order to use the property of black edges map[u] ^ map[v] = 1

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
ll multiplier = 1; // for computing the hash
int possible_pairs = 0, valid_pairs = 0; // storing the number of valid bibubble pairs

vector<int> stack_trace; // for storing the vertices in the stack during dfs traversal
vector<int> rm_cnt; // vector storing the number of brackets ending at that node
vector<pii> canonical_sese; // for storing the canonical sese pairs
unordered_map<ll, int> st; // (edge, size) -> hashed into a ll key
vector<bracketlist*> bl; // vector storing the bracket list for different nodes
vector<vector<edge*>> remove_brackets;

ll get_key(int id, int size){
    return size * multiplier + id;
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
    printArgs("sese:", u, parent);
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
        // eid.pb({min(w, u), max(w, u)});
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
                printArgs("cycle equivalent pair:", u, w);
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
vector<int> bb_comp; // bibubble-component to which a node belongs to
vector<vector<int>> bb_nodes; // storing the nodes for different bibubbles
// vector<bool> brute_valid; // checking whether a visited node belongs to a valid bibubble during the brute check
vector<bool> valid; // for marking whether a potential bi-bubble is valid
// vector<short> brute_visit; // checking whether a node is visit during the brute check for valid bibubbles
vector<int> aux_cc_comp; // connected-component in the auxiliary graph to which a node belongs to
vector<int> order; // will be a sorted list of G's vertices by exit time
vector<vector<int>> ag; // directed graph
vector<vector<int>> ag_rev; // create adjacency list of G^T   
queue<int> q; // queue for checking the found bibubbles using brute force

string get_label(int x, int y){
    x = id_in_original_graph[x]; y = id_in_original_graph[y];
    return (ilmap[x >> 1] + " " + (x & 1 ? "+" : "-") + " " + ilmap[y >> 1] + " " + (y & 1 ? "-" : "+"));
}

bool end_gene(int& u, pii& rs){
    return u == rs.F || u == dual[rs.F] || u == rs.S || u == dual[rs.S];
}

void add_directed_edges_from_gene_endpoints(int& id1, int& id2){
    ag[dual[id1]].pb(id2); ag[dual[id2]].pb(id1);
}

void upd_node(int u, int id){
    mark[u] = true;
    bb_nodes[id].pb(u);
    bb_comp[u] = id;
}

void mark_bb_nodes(int u, pii& rs, int id){
    // if((u == rs.S || u == dual[rs.S]) && bb_nodes[id].size() != 0)return; // u == dual[rs.S] -> can go to y complement but then the pair is not cycle equivalent, if size = 0 implies hairpin loop
    if(u != rs.F){// only the directed genes \in \tilde{U} set
        upd_node(u, id);
    }

    for(pii child : g_processed[u]){
        int v = child.F;
        if(child.S >= tip_start[biedged_connected_comp[rs.S]] || end_gene(v, rs))continue; // don't consider the edges added because of tips
        if(mark[v]){
            if(bb_comp[v] != bb_comp[u])bb_nodes[id].pb(v); // if the two bibubbles are adjacent, then too start and end points of the bibubble will be added, however the complete set need not be added
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

// bool brute(pii rs, int type){
//     q.push(type == 0 ? rs.F : rs.S);

//     while(!q.empty()){
//         int u = q.front(); q.pop();
//         for(pii child : g_processed[u]){
//             int v = child.F;
//             if(end_gene(v, rs))continue;
//             if(brute_valid[bb_comp[v]]){
//                 if(type == 0){
//                     if(brute_visit[v] == 0){
//                         q.push(dual[v]); brute_visit[v] = 1; brute_visit[dual[v]] = 1;
//                     }
//                 }else{
//                     if(brute_visit[v] == 0)return false;
//                     else if(brute_visit[v] == 1){
//                         q.push(dual[v]); brute_visit[v] = 2; brute_visit[dual[v]] = 2;
//                     }
//                 }
//             }else return false;
//         }
//     }
//     return true;
// }
// ****************************************************************************************************************************************************** 

int main(int argc, char* argv[])
{   
    // ************************************
    // *** IO + data preparation ***
    // ************************************
    {
        ios_base::sync_with_stdio(false); cin.tie(0); cout.tie(0); // Fast IO
        if(argc != 3){
            return 1;
        }
        inputpath = argv[1];
        outputdir = argv[2];
        if(!fs::exists(outputdir))fs::create_directories(outputdir);
        
        summarypath = outputdir + "/input_summary_debug.txt";
        freopen(summarypath.c_str(), "w", stdout);

        get_ne();

        cout << "Number of nodes in the input: " << n << endl;
        cout << "Number of edges in the input: " << edges << endl;

        cout << "HI\n";
        return 0;
        
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
                            if(cnt_vu == 0 && !has_self_loop[v]){ // tip case
                                type_edge[i] = '0';
                            }else{
                                n_processed++;
                                type_edge[i] = '1';
                            }
                        }else if(splitfrom_v){
                            if(cnt_uv == 0 && !has_self_loop[u]){
                                type_edge[i] = '0';
                            }else{
                                n_processed++;
                                type_edge[i] = '2';
                            }                        
                        }else{
                            type_edge[i] = '4';
                        }
                    }else type_edge[i] = '0';
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

                // printGraph(g_processed);

                // ** Clearing the memory **
                g.clear(); type_edge.clear();
                // vector<vector<pii>>().swap(g); vector<char>().swap(type_edge);
            }
        }
        
        // ************************************
        // *** Finding number of connected components in the processed graph***    
        // count of vertices = n_processed
        // count of edges = n (no of genes)
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
                if(tips[it].size() <= 1)tip_start[it] = __LONG_LONG_MAX__; 
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
        
        // printEdges();

        // ************************************
        // *** SESE ***
        // ************************************
        {
            remove_brackets.resize(n_processed);
            bl.clear();
            fill(mark.begin(), mark.end(), false);
            for(int it = 0; it < after_initialisation_comp_cnt; it++){
                multiplier = 1;
                while(multiplier <= backedge_cnt[it]){
                    multiplier *= 10;
                }
                st.clear();
                sese(Start[it], -1); // graph g is not changed
                // assert(stack_trace.size() == 0); 

                // corner case
                if(tips[it].size() == 1){
                    canonical_sese.pb({dual[tips[it][0]], dual[tips[it][0]]});
                }
            }

            // ** Clearing the memory **
            bl.clear(); remove_brackets.clear(); Start.clear(); depth.clear(); 
            backedge_cnt.clear(); tips.clear();
            // vector<bracketlist*>().swap(bl); vector<vector<edge*>>().swap(remove_brackets); vector<int>().swap(Start); vector<int>().swap(depth);  
            // vector<int>().swap(backedge_cnt);  vector<vector<int>>().swap(tips);
        }
        
        // ************************************
        // *** Printing result ***
        // ************************************
        {
            summarypath = outputdir + "/possible_bibubble_debug.txt";
            freopen(summarypath.c_str(), "w", stdout);

            possible_pairs = canonical_sese.size();
            cout << "Total possible canonical cycle equivalent pairs found: " << possible_pairs << endl;
            
            for(int i = 0; i < possible_pairs; i++){
                pii rs = canonical_sese[i];
                cout << get_label(rs.F, rs.S) << endl;                  
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
                
                for(int i = 0; i < possible_pairs; i++){ // space requirement is linear since the inner-most nested bubbles will appear on the top
                    pii rs = canonical_sese[i];
                    // upd_node(id_in_original_graph[rs.F] ^ 1, i); // will create issues when two bibubbles are adjacent
                    mark_bb_nodes(rs.F, rs, i); // do not add the end points, checking only in one direction
                }

                // for(int i = 0; i < possible_pairs; i++){
                //     printArgs("bb_nodes", i, ":");
                //     printVector(bb_nodes[i]);
                // }
            }
        }

        valid.resize(possible_pairs, true);
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
                    if(!((aux_cc_comp[u] == aux_cc_comp[rs.F] || aux_cc_comp[u] == aux_cc_comp[dual[rs.F]]) && valid[bb_comp[u]])){
                        valid[i] = false;
                        valid_pairs--;
                        break;
                    }
                // }
            }
            if(bb_nodes[i].size() == 0){
                valid[i] = false;
                valid_pairs--;
            }
        }

        // ************************************
        // *** Printing the result ***
        // ************************************        
        {
            // ** Clearing the memory **
            tip_start.clear(); 
            // vector<ll>().swap(tip_start);

            summarypath = outputdir + "/valid_bibubble_GO_debug.txt";
            freopen(summarypath.c_str(), "w", stdout);

            cout << "Total valid canonical cycle equivalent pairs found: " << valid_pairs << endl;
            for(int i = 0; i < possible_pairs; i++){
                if(!valid[i])continue;
                pii rs = canonical_sese[i];
                cout << get_label(rs.F, rs.S) << endl;
            }
        }
    }
} 