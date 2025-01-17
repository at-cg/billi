// Remove eid, scnt, free unused memory, implementation for different components

#include<iostream>
#include<fstream>
#include<filesystem>
#include<vector>
#include<map>
#include<set>
#include<unordered_map>
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

int n; // no of nodes (genes)
int comp_cnt = 0; // no of connected components in the input graph
int nodes = 0; // no of nodes in a given componenet in the input graph
ll multiplier = 1; // for computing the hash
int maxd; // max depth possible = total number of nodes
int possible_pairs = 0, valid_pairs = 0; // storing the number of valid bibubble pairs
// int End = -1; // for connecting the tips
int tot_grey = 0; // for storing total number of grey edges, will help in accessing the brackets

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
    bracketlist():sz(0), d1(maxd), d2(maxd + 10), start(nullptr), end(nullptr){}
    bracketlist(int lsz, int depth1, int depth2, edge* pstart, edge* pend): sz(lsz), d1(depth1), d2(depth2), start(pstart), end(pend){}
};

vector<bool> mark; // for marking back edges
vector<bool> valid; // for marking whether a potential bi-bubble is valid
vector<int> stack_trace; // for storing the vertices in the stack during dfs traversal
vector<int> depth; // for storing the depth of the vertices in the spanning tree
vector<int> first_node; // for storing a node in every connected-component in the biedged graph
vector<int> bi_cc_comp; // connected-component in the biedged graph to which a node belongs to 
vector<int> aux_cc_comp; // connected-component in the auxiliary graph to which a node belongs to
vector<int> bb_comp; // bibubble-component to which a node belongs to
vector<int> order; // will be a sorted list of G's vertices by exit time
vector<int> backedge_cnt; // upper bound on number of back edges
vector<int> Start; // starting point of search for each component in original graph
vector<ll> tip_start; // will store the id from which extra edges because of tips are added for each bi_cc_comp
vector<string> ilmap; // for storing the gene for a particular label
vector<pii> eid; // stores the edges
vector<pii> canonical_sese; // for storing the canonical sese pairs
vector<vector<int>> bb_nodes; // storing the nodes for different bibubbles
vector<vector<int>> ag; // directed graph
vector<vector<int>> ag_rev; // create adjacency list of G^T   
vector<vector<int>> tips; // vector for storing the vertices that represent the tip
vector<vector<pii>> g; // (node_id, grey_edge_id)
vector<vector<edge*>> remove_brackets;
unordered_map<ll, int> st; // (edge, size) -> hashed into a ll key
map<string, int> lmap; // for storing the label for a particular gene

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
        string s2 =  !(rs.S & 1) ? "-" : "+"; // signs are reversed
        add_directed_edges(id1, id2, s1, s2);
    }
}

void merge(bracketlist* bl1, bracketlist* bl2){
    if(!bl1->start){
        bl1->start = bl2->start;
        bl1->d1 = bl2->d1;
        bl1->d2 = bl2->d2;
    }else{
        bl1->end->front = bl2->start;
        bl2->start->back = bl1->end;
        if(bl1->d1 < bl2->d1){
            bl1->d2 = min(bl1->d2, bl2->d1);
        }else if(bl1->d1 > bl2->d1){
            bl1->d2 = min(bl1->d1, bl2->d2);
            bl1->d1 = bl2->d1;
        }else{
            bl1->d2 = min(bl1->d2, bl2->d2);
        }
        // else{
        //     bl1->d2 = min(bl1->d1, bl2->d2);
        //     bl1->d1 = bl2->d1;
        // }
    }
    bl1->end = bl2->end;
    bl1->sz += bl2->sz;
    delete(bl2);
}

void dfs_comp(int u){
    mark[u] = true;
    bi_cc_comp[u] = comp_cnt;
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

int find_unique(int u){
    set<int> s;
    for(pii child : g[u]){
        if((child.F^u) == 1)continue;
        s.insert(child.F);
        if(s.size() == 2)break;
    }
    return s.size();
}

bracketlist* sese(int u, int parent){
    stack_trace.pb(u);
    mark[u] = true;

    // printArgs(u, parent); 
    bracketlist* bl = new bracketlist(); 
    int cnt_back = 0;
    for(pii child : g[u]){
        int v = child.F;
        if(v == parent)continue;
        bracketlist* bl1;
        if(mark[v]){// back-edge -> will come first in the dfs traversal for the node at greater depth
            if(depth[v] > depth[u])continue; // front-edge
            assert(child.S != -1); // black edges can't be backedges
            edge* ed = new edge(child.S); // tot_grey count starts from 0
            bl1 = new bracketlist(1, depth[v], maxd, ed, ed);
            remove_brackets[v].pb(ed);
        }else bl1 = sese(v, u);
        if((u ^ v) == 1){// check for canonical sese if it is a black edge, won't be true for back edge
            // printArgs("Bracket list size on reaching the solid edge:", bl1->sz);
            ll key;
            if(bl1->sz > 0){
                int br = bl1->end->id;
                assert(br != -1); // not a black edge
                key = get_key(br, bl1->sz);
            }else{
                key = 0;
            }
            // if(u == 3996){
            //     scnt += 1;
            //     printArgs("HI", u, scnt, depth[u]);
            //     printBracketList(bl1);
            // }
            // if(v == 1301){
            //     scnt += 1;
            //     printArgs("HI", v, scnt, depth[v]);
            //     printBracketList(bl1);
            // }

            if(st.find(key) != st.end()){
                // printArgs("Found the canonical pair");
                // if(g[v].size() > 2 || g[st[key]].size() > 2){ // to avoid linear chains <- fails in case of multiple edges
                int w = st[key];
                if(find_unique(v) == 2 || find_unique(w) == 2 || (find_unique(v) == 1 && g[v][1].F != w) || (find_unique(w) == 1 && g[w][1].F != v)){ // to avoid linear chains
                    // printArgs("Found the contributing canonical pair:", v, st[key]);
                    canonical_sese.pb({v, w});
                }
            }
            st[key] = u; // will happen regardless you found something or not
        }
        if(bl1->start){
            // printArgs(u, "child", v);
            // printBracketList(bl); printBracketList(bl1);
            if(bl1->d1 < depth[u] && depth[v] > depth[u])cnt_back++; // an edge from u to its ancestor won't contribute to capping backedge
            // printArgs(bl1->d1, depth[u], depth[v], "count of back edges:", cnt_back);
            merge(bl, bl1); // merging brackets from the child nodes <- should come in the end since bl1 will be deleted after this operation
            // printBracketList(bl);
        }
    }

    // removing brackets from respective lists after computing the value for the outgoing edges from that node
    for(edge* ed : remove_brackets[u]){
        // printArgs("remove", eid[ed->id].F, eid[ed->id].S);
        // int min_depth = min(depth[eid[ed->id].F], depth[eid[ed->id].S]);
        // assert(min_depth == depth[u]);
        if(bl->d1 == depth[u])bl->d1 = maxd;
        else if(bl->d2 == depth[u])bl->d2 = maxd;
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
        // printBracketList(bl);
    }
   
    // only left with the brackets that go up the node u
    // adding capping backedge
    // edge to be added u -> stack_trace[bl->d2]
    if(cnt_back >= 2 && bl->d2 < maxd){// at least two childs should have backedges with depth lower than node u in order to add a capping backedge
        int w = stack_trace[bl->d2]; // stack_trace is 0-indexed and so are bl->d1 and bl->d2
        // edge need not be added to graph g
        // check if the edge already exists
        bool exists = false;
        for(pii child : g[u]){
            if(child.F == w){
                exists = true; break;
            }
        }
        // printArgs("Capping back-edge", u, w);
        if(!exists){
            // printArgs("Capping back-edge does not exist");
            edge* ed = new edge(tot_grey); tot_grey++;
            eid.pb({min(w, u), max(w, u)});
            bracketlist* bl1 = new bracketlist(1, depth[w], maxd, ed, ed);
            remove_brackets[w].pb(ed);
            merge(bl, bl1);
        }
    }

    // printArgs(u, parent, stack_trace.size(), bl->sz);
    stack_trace.rb();
    return bl;
}

void upd_node(int u, int id){
    mark[u] = true;
    bb_nodes[id].pb(u);
    bb_comp[u] = id;
}

void mark_bb_nodes(int u, int end, int id){
    // printArgs("u, end, id:", u, end, id);
    upd_node(u, id);
    if(u == end){
        upd_node(u^1, id);
    }
    for(pii child : g[u]){
        int v = child.F;
        if(child.S >= tip_start[bi_cc_comp[end]])continue; // don't consider the edges added because of tips
        // printArgs("u, v:", u, v, mark[u], mark[v]);
        if(mark[v]){
            if(bb_comp[v] != bb_comp[u])bb_nodes[id].pb(v); // if the two bibubbles are adjacent, then too start and end points of the bibubble will be added, however the complete gene maynot be added
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
    string summarypath = outputdir + "/summary.txt";
    freopen(summarypath.c_str(), "w", stdout);

    get_n();
    maxd = 2 * n + 100;

    // printArgs(lmap["BENKGBHG_01759"], lmap["Rv1676"]);
    // deg.resize(2 * n);

    // for a node x (0-indexed) in pangene graph - two nodes 2 * x (tail of arrow), 2 * x + 1 (head of arrow) are created in the bi-edged graph 
    // using vectors will be good coz we will have to add capping backedges as well
    g.resize(2 * n); // // +delta is for S if required -> not needed
    // adding black edges first <- important since don't want black edge to appear as a back or front edge
    for(int i = 0; i < n; i++){
        int n1 = 2 * i, n2 = 2 * i + 1;
        g[n1].pb({n2, -1}); g[n2].pb({n1, -1});
    }
    make_graph(); // adding grey edges

    // ************************************
    // Finding number of connected components    
    // ************************************
    mark.resize(2 * n);
    bi_cc_comp.resize(2 * n);

    for(int i = 0; i < 2 * n; i += 2){
        if(mark[i]){
            assert(mark[i + 1]);
            continue;
        }
        dfs_comp(i);
        first_node.pb(i);
        comp_cnt++;
    }
    // printArgs("Number of components:", comp_cnt);

    // Graph compression - removing unitigs --- Not implementing currently -> does not seem that useful

    // ************************************
    // Identifying tips
    // Basically the vertices with no of unique neighbors = 1
    // ************************************
    tips.resize(comp_cnt);
    for(int i = 0; i < 2 * n; i++){
        // if(find_unique(i) == 0){// only connected via a black edge
        if(g[i].size() == 1){
            tips[bi_cc_comp[i]].pb(i);
            assert((g[i][0].F ^ i) == 1); // checking whether the node is connected via a black edge
        }
    }
    // for(int i = 0; i < comp_cnt; i++){
    //     printArgs("Tips", i);
    //     printVector(tips[i]);
    // }
    // random_shuffle(tips.begin(), tips.end());

    // ************************************
    // Finding depth and number of backedges for individual components
    // ************************************
    tip_start.resize(comp_cnt, -1);
    backedge_cnt.resize(comp_cnt);
    Start.resize(comp_cnt, -1);
    fill(mark.begin(), mark.end(), false);
    depth.resize(2 * n);
      
    for(int it = 0; it < comp_cnt; it++){   
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
        
        // Finding an upperbound on the number of backedges
        nodes = 0;
        dfs_depth(Start[it], -1, it);
        backedge_cnt[it] += nodes; // capping backedges
    }

    // for(int i = 0; i < eid.size(); i++){
    //     printArgs(i, eid[i].F, eid[i].S);
    // }
    
    // ************************************
        // *** Finding Possible Bibubble Pairs ***
    // ************************************
    // SESE
    // An important observation is that bibubble ends won't turn as backedges
        
    remove_brackets.resize(2 * n);
    fill(mark.begin(), mark.end(), false);
    for(int it = 0; it < comp_cnt; it++){
        multiplier = 1;
        while(multiplier <= backedge_cnt[it]){
            multiplier *= 10;
        }
        // printArgs("Multiplier:", multiplier);

        st.clear();
        sese(Start[it], -1); // graph g is not changed
        assert(stack_trace.size() == 0); 
    }
    
    possible_pairs = canonical_sese.size();
    cout << "Total possible canonical cycle equivalent pairs found: " << possible_pairs << endl;

    fill(mark.begin(), mark.end(), false);
    bb_comp.resize(2 * n, -1);
    bb_nodes.resize(possible_pairs);
    valid.resize(possible_pairs);

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
            for(int u : bb_nodes[i]){
                assert(bb_comp[u] != -1);
                // printArgs(rs.F, "-", aux_cc_comp[rs.F], rs.F^1, "-", aux_cc_comp[rs.F^1], u, "-", aux_cc_comp[u]);
                if(u != rs.F && u != (rs.F ^ 1) && u != rs.S && u != (rs.S ^ 1)){ // check only for internal nodes
                    if(!((aux_cc_comp[u] == aux_cc_comp[rs.F] || aux_cc_comp[u] == aux_cc_comp[rs.F^1]) && valid[bb_comp[u]])){
                        valid[i] = false;
                        valid_pairs--;
                        break;
                    }
                }
            }
        }
    }

    cout << "Total valid canonical cycle equivalent pairs found: " << valid_pairs << endl;
    for(int i = 0; i < possible_pairs; i++){
        if(!valid[i])continue;
        pii rs = canonical_sese[i];
        cout << get_label(rs.F, rs.S) << endl;
    }
} 
