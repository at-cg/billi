#include "commons.hpp"

// ****************************************************************************************************************************************************** 
// **** Generic ****
// ****************************************************************************************************************************************************** 
int tot_gray = 0; // for storing total number of gray edges, will help in accessing the brackets
vector<bool> mark; // for marking the vertices that have been visited
vector<string> ilmap; // for storing the gene for a particular label
vector<pii> edge_label; // edge vertex labels for the given edge id in the compacted graph
bool PRINT = false; // for debugging

string get_single_label(int& x, int ty){
    return ((ty == 0 ? (x & 1 ? ">" : "<") : (x & 1 ? "<" : ">")) + ilmap[x >> 1]);
}

string get_label(int x, int y){
    return (get_single_label(x, 0) + "\t" + get_single_label(y, 1)); 
}

void increase_stack() {
    const rlim_t kStackSize = RLIM_INFINITY;
    struct rlimit rl;

    if (getrlimit(RLIMIT_STACK, &rl) == 0) {
        rl.rlim_cur = kStackSize;
        if (setrlimit(RLIMIT_STACK, &rl) != 0) {
            std::cerr << "setrlimit failed\n";
        }
    }
}
// ****************************************************************************************************************************************************** 

// ****************************************************************************************************************************************************** 
// **** Helper functions ****
// ****************************************************************************************************************************************************** 
// void printBracketList(bracketlist* bl){
//     if(!PRINT)return;
//     // printArgs("d1:", bl->d1, "d2:", bl->d2);
//     cout << "Bracket List: ";
//     edge* it = bl->start;
//     while(it){
//         cout << it->id << " (" << get_single_label(edge_label[it->id].F, 0) << " " << get_single_label(edge_label[it->id].S, 0) << ") ";
//         it = it->front;
//     }cout << endl;
// }

void printVector(vector<int>& v, int ty){
    if(!PRINT)return;
    for(int x : v){
        cout << (ty == 0 ? to_string(x) : get_single_label(x, 0)) << " ";
    }
    cout << endl;
}

void printVector(vector<pii>& v, int ty){
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
int n, edges; // no of nodes (genes), no of edges
vector<bool> has_self_loop; // whether the node has a self loop
vector<pss> hap_walk; // for storing the haplotype walks
vector<vector<pii>> g; // (node_id, gray_edge_id)
map<string, int> lmap; // for storing the label for a particular gene

void get_ne(string inputpath) {
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

void get_walk(string inputpath) {
    ifstream f(inputpath);
    string line;

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
                if (tokens[0] == "W") {
                    hap_walk.pb(mp(tokens[1] + "#" + tokens[2], tokens[6]));
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

void make_graph(string inputpath) {
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
int cnt_black_edge = 0; // count of black edges in the compacted graph
vector<int> dual; // neighbor of a vertex connected by a black edge
queue<pii> q_compact; // queue used for clustering the edges in a linear chain

bool chk_vertex(int& id){
    return g[id].size() == 2 && g[id][0].F != g[id][1].F && !has_self_loop[id];
}

bool chk_for_compaction(int& id){
    if(chk_vertex(id)){
        int u = g[id][1].F;
        return chk_vertex(u);
    }
    return false;
}
// ****************************************************************************************************************************************************** 
