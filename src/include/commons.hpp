#pragma once

#include<iostream>
#include<fstream>
#include<filesystem>
#include<sys/resource.h>

#include<vector>
#include<string>
#include<map>
#include<set>
#include<unordered_map>
#include<algorithm>
#include<stack>
#include<queue>
#include<deque>
#include<regex>
#include<cassert>
#include<stdexcept>

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

// ****************************************************************************************************************************************************** 
// **** Generic ****
// ****************************************************************************************************************************************************** 
extern int maxd; // max depth possible
extern int tot_gray; // for storing total number of gray edges, will help in accessing the brackets
extern vector<bool> mark; // for marking the vertices that have been visited
extern vector<string> ilmap; // for storing the gene for a particular label
extern vector<pii> edge_label; // edge vertex labels for the given edge id in the compacted graph
extern bool PRINT; // for debugging

string get_single_label(int& x, int ty);

string get_label(int x, int y);

void increase_stack();
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
void printArgs(const Args&... args){
    if(!PRINT)return;
    ((cout << args << " "), ...) << endl; // Fold expression
}

void printBracketList(bracketlist* bl);

void printVector(vector<int>& v, int ty = 0);

void printVector(vector<pii>& v, int ty = 0);

void printVector(vector<char>& v);

void printGraph(vector<vector<int>>& g);

void printGraph(vector<vector<pii>>& g);
// ****************************************************************************************************************************************************** 

// ****************************************************************************************************************************************************** 
// **** Input ****
// ****************************************************************************************************************************************************** 
extern int n, edges; // no of nodes (genes), no of edges
extern vector<bool> has_self_loop; // whether the node has a self loop
extern vector<vector<pii>> g; // (node_id, gray_edge_id)
extern map<string, int> lmap; // for storing the label for a particular gene

void get_ne(string inputpath);

void add_edge(vector<vector<pii>>& g, int id1, int id2);

void make_graph(string inputpath);
// ****************************************************************************************************************************************************** 

// ****************************************************************************************************************************************************** 
// **** Compacted Graph ****
// ****************************************************************************************************************************************************** 
extern int cnt_gray_edge; // count of gray edges in the compacted graph
extern int cnt_black_edge; // count of black edges in the compacted graph
extern vector<int> dual; // neighbor of a vertex connected by a black edge
extern queue<pii> q_compact; // queue used for clustering the edges in a linear chain

bool chk_vertex(int& id);

bool chk_for_compaction(int& id);
// ****************************************************************************************************************************************************** 
