// Remove vector eid, scnt, free unused memory, has_self_loop, clear bl vector

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

using namespace std;
namespace fs = std::filesystem;

fstream f;
string inputpath, outputdir;

// **** Generic ****
vector<bool> mark; // for marking the vertices that are have been visited
// *************************

// **** Input ****
int n, edges; // no of nodes (genes), no of edges
int before_initialisation_comp_cnt = 0; // no of connected components in the input graph before initialisation
vector<vector<pii>> g; // (node_id, grey_edge_id)
// ************************* 

// **** Processed Graph ****
int after_initialisation_comp_cnt = 0; // no of connected components in the input graph after initialisation
int n_processed; // no of nodes in the processed graph (at most {n + 2 * edges})
map<int, int> id_in_original_graph; // used to find the id in original graph, in order to use the property of black edges map[u] ^ map[v] = 1
vector<vector<pii>> g_processed; // (node_id, grey_edge_id)
// ************************* 


int nodes = 0; // no of nodes in a given component in the input graph
ll multiplier = 1; // for computing the hash
int maxd; // max depth possible = total number of nodes
int possible_pairs = 0, valid_pairs = 0; // storing the number of valid bibubble pairs
// int End = -1; // for connecting the tips
int tot_grey = 0; // for storing total number of grey edges, will help in accessing the brackets

// **** Finding bridges ****
int timer = 0; // for finding bridges in the uninitialised graph
vector<int> tin, low;
vector<bool> mark_bridge;

void dfs_bridge(int u, int parent){
    mark[u] = true;
    tin[u] = low[u] = timer++;

    for(pii child : g[u]){
        int v = child.F;
        if(v == parent)continue; // a multiedge b/w 2 nodes can also result in a bridge
        if(mark[v])low[u] = min(low[u], tin[v]);
        else{
            dfs_bridge(v, u);
            low[u] = min(low[u], low[v]);
            if(low[v] > tin[u] && (u ^ v == 1))mark_bridge[u >> 1] = true; // strictly greater sign (>) -- loop case
        }
    }
}
// *************************

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

// vector<bool> has_self_loop; // not required since an extra edge between a node to its parent will not lead to a cycle

vector<bool> brute_valid; // checking whether a visited node belongs to a valid bibubble during the brute check
vector<bool> valid; // for marking whether a potential bi-bubble is valid
vector<short> brute_visit; // checking whether a node is visit during the brute check for valid bibubbles
vector<int> stack_trace; // for storing the vertices in the stack during dfs traversal
vector<int> depth; // for storing the depth of the vertices in the spanning tree
vector<int> first_node; // for storing a node in every connected-component in the biedged graph
vector<int> biedged_connected_comp; // connected-component in the biedged graph to which a node belongs to 
vector<int> aux_cc_comp; // connected-component in the auxiliary graph to which a node belongs to
vector<int> bb_comp; // bibubble-component to which a node belongs to
vector<int> order; // will be a sorted list of G's vertices by exit time
vector<int> backedge_cnt; // upper bound on number of back edges
vector<int> Start; // starting point of search for each component in original graph
vector<ll> tip_start; // will store the id from which extra edges because of tips are added for each biedged_connected_comp
vector<string> ilmap; // for storing the gene for a particular label
vector<pii> eid; // stores the edges
vector<pii> canonical_sese; // for storing the canonical sese pairs
vector<vector<int>> bb_nodes; // storing the nodes for different bibubbles
vector<vector<int>> ag; // directed graph
vector<vector<int>> ag_rev; // create adjacency list of G^T   
vector<vector<int>> tips; // vector for storing the vertices that represent the tip
vector<vector<edge*>> remove_brackets;
unordered_map<ll, int> st; // (edge, size) -> hashed into a ll key
map<string, int> lmap; // for storing the label for a particular gene
vector<bracketlist*> bl; // vector storing the bracket list for different nodes
queue<int> q; // queue for checking the found bibubbles using brute force

template <typename... Args>
void printArgs(Args... args){
    // return;
    ((cout << args << " "), ...) << endl; // Fold expression
}

void printBracketList(bracketlist* bl){
    // return;
    // printArgs("d1:", bl->d1, "d2:", bl->d2);
    cout << "Bracket List: ";
    edge* it = bl->start;
    while(it){
        cout << it->id << " ";
        it = it->front;
    }cout << endl;
}

void printVector(vector<int>& v){
    for(int x : v){
        cout << x << " ";
    }
    cout << endl;
}

void printGraph(vector<vector<int>>& g){
    int sz = g.size();
    for(int i = 0; i < sz; i++){
        cout << i << ": ";
        printVector(g[i]);
    }
}

ll get_key(int id, int size){
    return size * multiplier + id;
}

inline string get_label(int x, int y){
    return (ilmap[x >> 1] + " " + (x & 1 ? "+" : "-") + " " + ilmap[y >> 1] + " " + (y & 1 ? "-" : "+"));
}

void get_n(){
    f.open(inputpath, ios::in);
    string line;
    regex strip("^\\s+|\\s+$"), split("\\t");
    vector<string> tokens;
    n = 0; 
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
                lmap[tokens[1]] = n; // 0-indexed
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

void get_e(){
    f.open(inputpath, ios::in);
    string line;
    regex strip("^\\s+|\\s+$"), split("\\t");
    vector<string> tokens;
    edges = 0; 
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
                edges++;
            }
        }
        f.close();
    }
}

void add_directed_edges(int& id1, int& id2, string s1, string s2){
    // printArgs(id1, id2, s1, s2);
    if(s1 == s2){ // >x denotes +x 
        if(s1 == "+"){
            ag[id1].pb(id2); ag[id2^1].pb(id1^1);
        }else{
            ag[id2].pb(id1); ag[id1^1].pb(id2^1);
        }
    }else{
        if(s1 == "+"){
            ag[id1].pb(id2^1); ag[id2].pb(id1^1);
        }else{
            ag[id1^1].pb(id2); ag[id2^1].pb(id1);
        }
    }
}

void add_edge(int& id1, int& id2){
    g[id1].pb({id2, tot_grey}); g[id2].pb({id1, tot_grey}); // 0-indexed
    eid.pb({min(id1, id2), max(id1, id2)});
    tot_grey++;
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
                add_edge(id1, id2);
                // printArgs(id1, id2, tot_grey);

                // if(n1 == n2 && s1 == s2){
                //     has_self_loop[id1] = has_self_loop[id2] = true;
                // }
            }
        }
        f.close();
    }
}

void make_auxillary_graph(){
    f.open(inputpath, ios::in);
    string line;
    regex strip("^\\s+|\\s+$"), split("\\t");
    vector<string> tokens;

    // using the input edges
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
                // if(n1 == n2)continue; // can be skipped here since self loops can help turn back through the same node
                
                // making sure edges are added only in one direction
                if(n1 > n2){
                    swap(n1, n2); 
                    if(s1 == "+")s1 = "-";
                    else s1 = "+";
                    if(s2 == "+")s2 = "-";
                    else s2 = "+";
                    swap(s1, s2);
                }
                
                // n1 and n2 are 0-indexed
                int id1 = n1 << 1, id2 = n2 << 1;
                add_directed_edges(id1, id2, s1, s2);
            }
        }
        f.close();
    }

    // using the extra added edges
    for(pii rs : canonical_sese){
        int id1 = (rs.F | 1) ^ 1, id2 = (rs.S | 1) ^ 1;
        string s1 = (rs.F & 1) ? "-" : "+";
        string s2 = (rs.S & 1) ? "+" : "-"; // signs are reversed
        add_directed_edges(id1, id2, s1, s2);
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
    // printArgs("y3");
    // delete(bl2);
    bl2->merged = true;
    // printArgs("y4");
}

void dfs_comp(int u){ // Finding connected components in the biedged graph
    mark[u] = true;
    biedged_connected_comp[u] = before_initialisation_comp_cnt;
    for(pii child : g[u]){
        int v = child.F;
        if(mark[v])continue;
        dfs_comp(v);
    }
}

void dfs_depth(int u, int parent, int comp){
    if(parent != -1 && !mark[u]){
        depth[u] = depth[parent] + 1;
        // printArgs(u, ":", depth[u]);
    }
    if(!mark[u])nodes++;
    mark[u] = true;
    for(pii child : g[u]){
        int v = child.F;
        if(mark[v]){
            if(depth[v] < depth[u] && v != parent)backedge_cnt[comp]++;
            continue;
        }
        dfs_depth(v, u, comp);
    }
}

int find_unique_excluding_selfloop(int u, int v){
    set<int> s;
    for(pii child : g[u]){
        if((child.F^u) <= 1 || child.F == v)continue;
        s.insert(child.F);
        break;
    }
    return s.size();
}

int find_unique_including_selfloop(int u){
    set<int> s;
    for(pii child : g[u]){
        if((child.F^u) == 1)continue;
        s.insert(child.F);
        if(s.size() == 2)break;
    }
    return s.size();
}

void sese(int u, int parent){
    // printArgs(u, parent, "in");

    stack_trace.pb(u);
    mark[u] = true;

    int h0 = maxd, h1 = maxd, h2 = maxd;

    for(pii child : g[u]){
        int v = child.F;
        if(v == parent || v == u)continue; // v = u won't contribute anything and v == parent won't contribute to backedge
        if(mark[v]){// back-edge -> will come first in the dfs traversal for the node at greater depth
            if(depth[v] > depth[u])continue; // front-edge -- multi edges will be taken care of here
            assert(child.S != -1); // black edges can't be backedges
            h0 = min(h0, depth[v]);
        }else sese(v, u);
    }

    // printArgs(u, parent, "out");

    bl[u] = new bracketlist(); 

    // printArgs(u, parent, "h01");
    for(pii child : g[u]){
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
        // printArgs("y1", u, v);
        if(!bl[v] -> merged){
            // printBracketList(bl[u]); printBracketList(bl[v]);
            merge(bl[u], bl[v]);
            // printBracketList(bl[u]);
        }
        // printArgs("y2");
    }
    // printArgs(u, parent, "h02");
    bl[u]->d = min(h0, h1);

    // printArgs(h0, h1, h2);

    // removing brackets from respective lists
    // printArgs(u, parent, "h11");
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
    // printArgs(u, parent, "h12");

    // pushing back edges from node u
    // printArgs(u, parent, "h21");
    for(pii child : g[u]){
        int v = child.F;
        if(depth[v] >= depth[u] - 1)continue;
        edge* ed = new edge(child.S); // tot_grey count starts from 0
        remove_brackets[v].pb(ed);
        merge(bl[u], new bracketlist(1, depth[v], ed, ed));
    }
    // printArgs(u, parent, "h22");

    // capping back edge
    // printArgs(u, parent, "h31");
    if(h2 < h0){
        // printArgs("finding capping backedge");
        int w = stack_trace[h2];
        edge* ed = new edge(tot_grey); tot_grey++; // need not modify g
        eid.pb({min(w, u), max(w, u)});
        remove_brackets[w].pb(ed);
        merge(bl[u], new bracketlist(1, depth[w], ed, ed));
        // printArgs("capping back edge:", u, w);
    }
    // printArgs(u, parent, "h32");

    // finding cycle equivalence
    // printArgs(u, parent, "h41");
    if(parent != -1 && ((parent ^ u) == 1)){
        // printArgs("finding cycle equivalence:", u, parent);
        ll key;
        if(bl[u]->sz > 0){
            int br = bl[u]->end->id;
            assert(br != -1); // not a black edge
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
    // printArgs(u, parent, "h42");

    stack_trace.rb();
}

void upd_node(int u, int id){
    mark[u] = true;
    bb_nodes[id].pb(u);
    bb_comp[u] = id;
}

void mark_bb_nodes(int u, int end, int id){
    // if(id <= 5)printArgs("u, end, id:", u, end, id);

    upd_node(u, id);
    if(u == end){
        upd_node(u^1, id);
    }
    for(pii child : g[u]){
        int v = child.F;
        if(child.S >= tip_start[biedged_connected_comp[end]])continue; // don't consider the edges added because of tips
        // if(id <= 5)printArgs("u, v:", u, v, mark[u], mark[v]);
        if(mark[v]){
            if(bb_comp[v] != bb_comp[u])bb_nodes[id].pb(v); // if the two bibubbles are adjacent, then too start and end points of the bibubble will be added, however the complete set maynot be added
            continue;
        }
        mark_bb_nodes(v, end, id);
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
    for (int i = 0; i < 2 * n; i++)
        if (!mark[i])
            scc_dfs(i, ag, order);

    ag_rev.resize(2 * n);
    
    for (int v = 0; v < 2 * n; v++)
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
}

bool end_gene(int u, pii rs){
    return u == rs.F || u == (rs.F ^ 1) || u == rs.S || u == (rs.S ^ 1);
}

bool brute(pii rs, int type){
    q.push(rs.F);

    while(!q.empty()){
        int u = q.front(); q.pop();
        for(pii child : g[u]){
            int v = child.F;
            if(end_gene(v, rs))continue;
            if(brute_valid[bb_comp[v]]){
                if(type == 0){
                    if(brute_visit[v] == 0){
                        q.push(v^1); brute_visit[v] = 1; brute_visit[v^1] = 1;
                    }
                }else{
                    if(brute_visit[v] == 0)return false;
                    else if(brute_visit[v] == 1){
                        q.push(v^1); brute_visit[v] = 2; brute_visit[v^1] = 2;
                    }
                }
            }else return false;
        }
    }
    return true;
}

int main(int argc, char* argv[])
{   
    // ************************************
    // *** IO + data preparation ***
    // ************************************
    
    ios_base::sync_with_stdio(false); cin.tie(0); cout.tie(0); // Fast IO
    if(argc != 3){
        return 1;
    }
    inputpath = argv[1];
    outputdir = argv[2];
    if(!fs::exists(outputdir))fs::create_directories(outputdir);
    
    string summarypath = outputdir + "/input_summary.txt";
    freopen(summarypath.c_str(), "w", stdout);

    get_n();

    cout << "Number of vertices: " << n << endl;

    get_e();

    cout << "Number of edges: " << edges << endl;

    maxd = 2 * n + 100;

    // for a node x (0-indexed) in pangene graph - two nodes 2 * x (tail of arrow), 2 * x + 1 (head of arrow) are created in the bi-edged graph 
    // using vectors will be good coz we will have to add capping backedges as well
    g.resize(2 * n); // +delta is for S if required -> not needed (choose S as one of the tip ends)
    // has_self_loop.resize(2 * n);
    // adding black edges first <- important since don't want black edge to appear as a back or front edge
    for(int i = 0; i < n; i++){
        int n1 = 2 * i, n2 = 2 * i + 1;
        g[n1].pb({n2, -1}); g[n2].pb({n1, -1});
    }
    make_graph(); // adding grey edges

    // ************************************
    // *** Initialisation ***    
    // ************************************

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
                    assert(mark[i + 1]);
                    continue;
                }
                dfs_bridge(i, -1);
            }

            // ** Clearing the memory **
            tin.clear(); low.clear();
            vector<int>().swap(tin); vector<int>().swap(low);
        }

        // ************************************
        // *** Splitting the graph from the bridge point ***  // done carefully such that unitig compression is not required
        // ************************************
        {   
            // ** Finding new n **
            n_processed = n;

            for(int i = 0; i < n; i++){ // iterate over all the black edges
                if(mark_bridge[i]){
                    int u = 2 * i, v = u + 1;

                }
            }
            
            // ** Clearing the memory **
            mark_bridge.clear(); 
            vector<bool>().swap(mark_bridge);
        }
        
        // ************************************
        // *** Finding number of connected components ***    
        // ************************************
        {
            fill(mark.begin(), mark.end(), false);
            biedged_connected_comp.resize(2 * n);

            for(int i = 0; i < 2 * n; i += 2){
                if(mark[i]){
                    assert(mark[i + 1]);
                    continue;
                }
                dfs_comp(i);
                first_node.pb(i);
                before_initialisation_comp_cnt++;
            }
        }

    // ************************************
    // *** Identifying tips ***
    // ************************************
    tips.resize(before_initialisation_comp_cnt);
    for(int i = 0; i < 2 * n; i++){
        // only connected via a black edge
        if(find_unique_including_selfloop(i) == 0){ // much stricter condition than below
        // if(g[i].size() == 1){
            tips[biedged_connected_comp[i]].pb(i);
            // assert((g[i][0].F ^ i) == 1); // checking whether the node is connected via a black edge
        }
    }
    // for(int i = 0; i < before_initialisation_comp_cnt; i++){
    //     printArgs("Tips", i);
    //     printVector(tips[i]);
    // }
    // random_shuffle(tips.begin(), tips.end());

    // ************************************
    // *** Finding depth and number of backedges for individual components ***
    // ************************************
    tip_start.resize(before_initialisation_comp_cnt, -1);
    backedge_cnt.resize(before_initialisation_comp_cnt);
    Start.resize(before_initialisation_comp_cnt, -1);
    fill(mark.begin(), mark.end(), false);
    depth.resize(2 * n);
    
    cout << "Tips size: ";
    for(int it = 0; it < before_initialisation_comp_cnt; it++){   
        cout << tips[it].size() << " ";
        if(tips[it].size() == 0){ // there can be multiple components here, need to take care of in the dfs 
            Start[it] = first_node[it]; // no extra edges required
            tip_start[it] = __LONG_LONG_MAX__;
        }else{
            Start[it] = tips[it][0]; //End = n; n++;
            // if(tips.size() > 1){// else no need to add an extra edge
                for(int i = 1; i < tips[it].size(); i++){
                    if(tip_start[it] == -1)tip_start[it] = tot_grey; // marking the starting id for the edges added because of tips
                    add_edge(Start[it], tips[it][i]);
                    // add_edge(End, tips[it][i]);
                }
                // add_edge(Start, End);
            // }
        }
        
        // Finding an upperbound on the number of backedges <- only these will contribute to the brackets
        nodes = 0;
        dfs_depth(Start[it], -1, it);
        backedge_cnt[it] += nodes; // capping backedges
    }

    cout << endl;

    // for(int i = 0; i < eid.size(); i++){
    //     printArgs(i, eid[i].F, eid[i].S);
    // }
    
    // ************************************
    // *** Finding Possible Bibubble Pairs ***
    // ************************************
    // SESE
    // An important observation is that bibubble ends won't turn as backedges
        
    remove_brackets.resize(2 * n);
    bl.resize(2 * n);
    fill(mark.begin(), mark.end(), false);
    for(int it = 0; it < before_initialisation_comp_cnt; it++){
        multiplier = 1;
        while(multiplier <= backedge_cnt[it]){
            multiplier *= 10;
        }
        // printArgs("Multiplier:", multiplier);

        st.clear();
        sese(Start[it], -1); // graph g is not changed
        assert(stack_trace.size() == 0); 
    }
    
    // return 0;

    possible_pairs = canonical_sese.size();
    cout << "Total possible canonical cycle equivalent pairs found: " << possible_pairs << endl;

    fill(mark.begin(), mark.end(), false);
    bb_comp.resize(2 * n, -1);
    bb_nodes.resize(possible_pairs);
   
    for(int i = 0; i < possible_pairs; i++){ // space requirement is linear since the inner-most nested bubbles will appear on the top
        pii rs = canonical_sese[i];
        cout << get_label(rs.F, rs.S) << endl;
        upd_node(rs.F^1, i);
        mark_bb_nodes(rs.F, rs.S, i);
        // printVector(bb_nodes[i]);
    }

    // ************************************
    // *** Finding Valid Bibubble Pairs ***
    // ************************************
    
    valid.resize(possible_pairs);
    if(possible_pairs != 0){
        fill(valid.begin(), valid.end(), true);
        aux_cc_comp.resize(2 * n, -1);
        ag.resize(2 * n);
        valid_pairs = possible_pairs;
        // need not connect tips here
        make_auxillary_graph(); // read the file again since this time the edges will be directed
        // printGraph(ag);
        scc(); // ref - https://cp-algorithms.com/graph/strongly-connected-components.html extended to multigraph
        for(int i = 0; i < possible_pairs; i++){
            pii rs = canonical_sese[i];
            // for(pii child : g[rs.F]){
            //     if(child.F == rs.S){ // removing adjacent pairs that might have been added because of edge redundancy
            //         valid[i] = false;
            //         valid_pairs--;
            //         break;
            //     }
            // }
            // if(!valid[i])continue;
            // if(i <= 5)printArgs(">", rs.F, rs.S);
            for(int u : bb_nodes[i]){
                // if(i <= 5){
                //     cout << u << " ";
                // }
                assert(bb_comp[u] != -1);
                // printArgs(rs.F, "-", aux_cc_comp[rs.F], rs.F^1, "-", aux_cc_comp[rs.F^1], u, "-", aux_cc_comp[u]);
                if(!end_gene(u, rs)){ // check only for internal nodes
                    if(!((aux_cc_comp[u] == aux_cc_comp[rs.F] || aux_cc_comp[u] == aux_cc_comp[rs.F^1]) && valid[bb_comp[u]])){
                        valid[i] = false;
                        valid_pairs--;
                        break;
                    }
                }
            }
            // if(i <= 5)cout << endl;
        }
    }

    cout << "Total valid canonical cycle equivalent pairs found: " << valid_pairs << endl;
    for(int i = 0; i < possible_pairs; i++){
        if(!valid[i])continue;
        pii rs = canonical_sese[i];
        cout << get_label(rs.F, rs.S) << endl;
    }

    // ************************************
    // *** Checking for FP using brute  ***
    // ************************************
    brute_visit.resize(2 * n);
    brute_valid.resize(possible_pairs, true);
    
    cout << "FP bibubbles found: " << endl;
    for(int i = 0; i < possible_pairs; i++){
        if(!valid[i])continue;
        while(!q.empty())q.pop();
        pii rs = canonical_sese[i];

        brute_valid[i] = brute(rs, 0);
        if(brute_valid[i])brute_valid[i] = brute_valid[i] & brute(rs, 1);
        if(!brute_valid[i]){
            cout << get_label(rs.F, rs.S) << endl;
        }
    }
} 

// todo - chk for vector initialisations