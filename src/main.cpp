#include<bits/stdc++.h>
#include<filesystem>
#define ll long long int
#define pb push_back
#define rb pop_back
#define ti tuple<int, int, int>
#define pii pair<int, int>
#define piii pair<int, pii>
#define pli pair<ll, int>
#define pll pair<ll, ll>
#define mp make_pair
#define mt make_tuple
#define F first
#define S second

using namespace std;
namespace fs = std::filesystem;

fstream f;
string inputpath, outputdir;
// vector<int> deg;
vector<vector<int>> g;
vector<bool> visit; // for marking back edges

int Ss = -1, Se = -1; // for connecting the tips

// Not required since the node labels are int's only in the gfa file 
// map<string, int> lmap; // for storing the label for a particular gene
// vector<int> ilmap; // for storing the gene for a particular label

struct node{
    int low, high; // a bracket can uniquely be identified by the lower and higher vertices it connects to : low -> lower height
    node* front;
    node* back; // for doubly linked list implementation -> helps in O(1) deletion

    node(int l, int h):low(l), high(h), front(nullptr), back(nullptr){}
};

inline int last_bit(string s){
    return s == "+" ? 0 : 1;
}

int get_n(){
    f.open(inputpath, ios::in);
    string line;
    regex strip("^\\s+|\\s+$"), split("\\t");
    vector<string> tokens;
    int n; // number of nodes in the pangene graph
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
                assert(tokens.size() == 3 || tokens.size() == 4);
                n = stoi(tokens[1]);
            }else if(tokens[0] == "L"){
                f.close();
                return n;
            }
        }
        f.close();
    }
    return -1;
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
                assert(tokens.size() == 6);
                int n1 = stoi(tokens[1]); string s1 = tokens[2]; 
                int n2 = stoi(tokens[3]); string s2 = tokens[4]; 
                n1--; n2--; // making 0-indexed
                int id1 = (n1 << 1) + last_bit(s1), id2 = (n2 << 1) + last_bit(s2);
                // deg[(n1 << 1) + last_bit(s1)]++; deg[(n2 << 1) + last_bit(s2)]++;
                g[id1].pb(id2); g[id2].pb(id1);
            }
        }
        f.close();
    }
}

vector<> sese(int u, int parent){
    visit[u] = true;
    for(int v : g[u]){
        if(v == parent)continue;
        if(visit[v]){

        }
    }
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

    int n = get_n();
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
    // Basically the vertices with deg = 1, also it's neighbor
    vector<int> tips; // vector for storing the vertices that represent the tip
    for(int i = 0; i < 2 * n; i++){
        if(g[i].size() == 1){
            tips.pb(i);
            assert(g[i][0] ^ i == 1); // checking whether the node is connected via a black edge
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
    visit.resize(2 * n);

    // Taking care of multiple components
    // for()
    sese(Ss, -1); 
} 
