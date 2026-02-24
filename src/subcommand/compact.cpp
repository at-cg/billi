#include "compact.hpp"

// ****************************************************************************************************************************************************** 
// **** Compacted Graph ****
// ****************************************************************************************************************************************************** 
vector<int> side; // for finding which side the edge belongs to
set<int> vertex_compacted; // stores the vertex id's in the compacted graph
set<int> edge_compacted; // stores the edge id's in the compacted graph
set<string> gray_edges; // for storing the gray edges in the compacted graph

void print_compacted_graph(string inputpath, string outputpath) {
    ifstream f(inputpath);
    string line;

    freopen(outputpath.c_str(), "w", stdout);

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

            if (!tokens.empty()) {
                if (tokens[0] == "S") {
                    if(edge_compacted.find(lmap[tokens[1]]) != edge_compacted.end()){
                        cout << line << endl;
                    }
                } else if (tokens[0] == "L") {
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

                    int min1 = min(id1, dual[id1]);
                    int min2 = min(id2, dual[id2]);

                    if(vertex_compacted.find(id1) != vertex_compacted.end() && vertex_compacted.find(id2) != vertex_compacted.end()){
                        tokens[1] = ilmap[min1 >> 1];
                        tokens[3] = ilmap[min2 >> 1];
                        
                        if(side[id1] == 1){
                            if(id1 == min1)tokens[2] = min1 & 1 ? "+" : "-";
                            else tokens[2] = min1 & 1 ? "-" : "+";
                        }else{
                            if(id1 == min1)tokens[2] = min1 & 1 ? "+" : "-";
                            else tokens[2] = min1 & 1 ? "+" : "-";
                        }

                        if(side[id2] == 1){
                            if(id2 == min2)tokens[4] = min2 & 1 ? "+" : "-";
                            else tokens[4] = min2 & 1 ? "+" : "-";
                        }else{
                            if(id2 == min2)tokens[4] = min2 & 1 ? "-" : "+";
                            else tokens[4] = min2 & 1 ? "+" : "-";
                        }

                        if(tokens[1] == tokens[3]){
                            if(tokens[2] == tokens[4] && tokens[2] == "-"){
                                tokens[2] = tokens[4] = "+";
                            }
                        }
                        string check = tokens[1] + tokens[2] + tokens[3] + tokens[4];
                        if(gray_edges.find(check) != gray_edges.end()){
                            continue;
                        }else{
                            gray_edges.insert(check);
                        }
                        
                        for(int i = 0; i < tokens.size(); i++){
                            cout << tokens[i] << (i == tokens.size() - 1 ? '\n' : '\t');
                        }
                    }
                } else if (tokens[0] == "W"){
                    string walk = tokens[6];
                    string compact_walk = "";
                    string id = "";
                    int node_id = -1;
                    for(char c : walk){
                        if(c == '>' || c == '<'){
                            if(id.length() > 0){
                                node_id = lmap[id.substr(1)];
                                if(edge_compacted.find(node_id) != edge_compacted.end()){
                                    compact_walk += id;       
                                }
                            }
                            id = c;
                        }else{
                            id += c;
                        }
                    }
                    if(id.length() > 0){
                        node_id = lmap[id.substr(1)];
                        if(edge_compacted.find(node_id) != edge_compacted.end()){
                            compact_walk += id;       
                        }
                    }

                    tokens[6] = compact_walk;
                    for(int i = 0; i < tokens.size(); i++){
                        cout << tokens[i] << (i == tokens.size() - 1 ? '\n' : '\t');
                    }
                }
            }
        }
        f.close();
    }
}
// ****************************************************************************************************************************************************** 

void run_compact(string inputpath, string outputpath)
{   
    // ************************************
    // *** IO + data preparation ***
    // ************************************
    {
        increase_stack(); // For setting the stack size to max_value
        ios_base::sync_with_stdio(false); cin.tie(0); cout.tie(0); // Fast IO

        get_ne(inputpath);

        cerr << "Number of vertices in the input bidirected graph: " << n << endl;
        cerr << "Number of edges in the input bidirected graph: " << edges << endl;

        // for a node x (0-indexed) in pangene graph - two nodes 2 * x (tail of arrow), 2 * x + 1 (head of arrow) are created in the bi-edged graph 
        // using vectors will be good coz we will have to add capping backedges as well
        g.resize(2 * n); // +delta is for S if required -> not needed (choose S as one of the tip ends)
        has_self_loop.resize(2 * n);
    
        // adding black edges first <- important since don't want black edge to appear as a back or front edge
        for(int i = 0; i < n; i++){
            int n1 = i << 1, n2 = n1 + 1;
            g[n1].pb({n2, -1}); g[n2].pb({n1, -1});
        }
        make_graph(inputpath); // adding gray edges
    }
    cerr << "Done reading input" << endl;
 
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
            mark.resize(2 * n, false);
            dual.resize(2 * n, -1);
            side.resize(2 * n, -1);

            for(int i = 0; i < 2 * n; i++){
                if(mark[i])continue;
   
                if(chk_for_compaction(i)){
                    q_compact.push({i ^ 1, 0}); q_compact.push({g[i][1].F ^ 1, 1});
                    int left = i, right = g[i][1].F;
                    mark[left] = mark[right] = true;
                    while(!q_compact.empty()){
                        pii pid = q_compact.front(); q_compact.pop();

                        if(chk_for_compaction(pid.F)){
                            mark[pid.F] = true;
                            int neigh = g[pid.F][1].F; mark[neigh] = true;
                            if(mark[neigh ^ 1])continue;
                            if(pid.S == 0){
                                left = neigh;
                                q_compact.push({left ^ 1, 0});
                            }else{
                                right = neigh;
                                q_compact.push({right ^ 1, 1});
                            }
                        }
                    }

                    int l = left ^ 1, r = right ^ 1;
                    dual[l] = r; dual[r] = l; 
                    side[l] = 0; side[r] = 1;
                    vertex_compacted.insert(l); vertex_compacted.insert(r);
                    edge_compacted.insert(min(l, r) >> 1);
                    cnt_black_edge++;

                    for(int j = 1; j < g[l].size(); j++){ // have removed parallel gray edges
                        int u = g[l][j].F;
                        if(mark[u] && u != r)continue; // u != r -> cyclic case
                        cnt_gray_edge++;
                    }
                    for(int j = 1; j < g[r].size(); j++){
                        int u = g[r][j].F;
                        if(mark[u])continue;
                        cnt_gray_edge++;
                    }
                    mark[l] = mark[r] = true; // done in the end so that self loops and parallel gray edges are not missed
                }
            }

            for(int i = 0; i < 2 * n; i++){
                if(mark[i])continue;

                if(!mark[i ^ 1]){
                    dual[i ^ 1] = i; dual[i] = i ^ 1;
                    side[i] = 0; side[i ^ 1] = 1;
                    vertex_compacted.insert(i); vertex_compacted.insert(i ^ 1);
                    edge_compacted.insert(i >> 1);
                    cnt_black_edge++;
                }
                for(int j = 1; j < g[i].size(); j++){
                    int u = g[i][j].F;
                    if(mark[u])continue;
                    cnt_gray_edge++;
                }
                mark[i] = true;
            }

            cerr << "Number of vertices in the bidirected graph after compaction: " << cnt_black_edge << endl;
            cerr << "Number of edges in the bidirected graph after compaction: " << cnt_gray_edge << endl;

            // ** Clearing the memory **
            g.clear();
        }
    }
    cerr << "Done graph cleaning" << endl;

    // ************************************
    // *** Printing output ***
    // ************************************
    {
        print_compacted_graph(inputpath, outputpath);
    }
    cerr << "Done printing compacted graph" << endl;
} 
