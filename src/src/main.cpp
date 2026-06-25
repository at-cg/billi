#include<CLI11.hpp>
#include<string>
#include<sstream>

#include "./subcommand/compact.hpp"
#include "./subcommand/decompose.hpp"

using namespace std;

string inputpath, outputpath;
bool print_reverse = false, use_exact = false, use_numeric = false, self_loop = false;

string get_cmdline(int argc, char** argv){
    ostringstream oss;
    for(int i = 0; i < argc; i++){
        if(i)oss << " ";
        oss << argv[i];
    }
    return oss.str();
}

int main(int argc, char* argv[])
{   
    string full_cmd = get_cmdline(argc, argv);

    CLI::App app{"Billi is a bubble-detection tool for pangenome graphs"};

    if(argc == 1){
        cout << app.help() << endl;
        return 0;    
    }

    auto compact = app.add_subcommand("compact", "Compact the input graph (in GFA format)");
    auto decompose = app.add_subcommand("decompose", "Decompose the input graph (in GFA format) into panbubbles");
    
    compact->add_option("-i, --input", inputpath, "Input file in GFA format")->required();
    compact->add_option("-o, --output", outputpath, "Output file in GFA format")->required();
    compact->add_flag("-r, --reverse", print_reverse, "Print edges in a reverse complement fashion as well (default: false)");
    compact->add_flag("-n, --numeric", use_numeric, "Use numeric id's as node labels (default: false)");
    compact->add_flag("-s, --self_loop", self_loop, "Retain self loops (default: false)");
    
    decompose->add_option("-i, --input", inputpath, "Input file in GFA format")->required();
    decompose->add_flag("-e, --exact", use_exact, "Use exact (slow) algorithm (default: heuristic)");

    app.require_subcommand(1);

    CLI11_PARSE(app, argc, argv);
    cerr << "Command line options:" << full_cmd << endl;

    if(*compact)run_compact(inputpath, outputpath, print_reverse, use_numeric, self_loop);
    else if(*decompose)run_decompose(inputpath, use_exact);
}
