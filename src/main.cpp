#include<CLI11.hpp>
#include<string>
#include<sstream>

#include "./subcommand/compact.hpp"
#include "./subcommand/decompose.hpp"

using namespace std;

string inputpath, outputpath;
bool print_reverse = false, use_exact = false, use_numeric = false, self_loop = false, iWalk = false;
int minAlleles = 0;

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

    auto compact = app.add_subcommand("compact", "Compute a compacted pangenome graph by collapsing non-branching linear paths");
    auto decompose = app.add_subcommand("decompose", "Report panbubbles and hairpins in the input pangenome graph");
    
    compact->add_option("-i, --input", inputpath, "Input file in GFA format")->required();
    compact->add_option("-o, --output", outputpath, "Output file in GFA format")->required();
    compact->add_flag("-r, --reverse", print_reverse, "Print edges in a reverse complement fashion as well (default: false)");
    compact->add_flag("-n, --numeric", use_numeric, "Modify node labels to numeric format (default: false)");
    compact->add_flag("-s, --self_loop", self_loop, "Retain self loops (default: false)");
    
    decompose->add_option("-i, --input", inputpath, "Input file in GFA format")->required();
    decompose->add_option("-m, --minAlleles", minAlleles, "Minimum number of alleles required for reporting a panbubble or hairpin (default: 0)");
    decompose->add_flag("-w, --iWalk", iWalk, "disable printing of alleles for faster execution (default: printed if GFA has W-lines)");
    decompose->add_flag("-e, --exact", use_exact, "Use exact algorithm instead of the fast heuristic approach (default: heuristic)");

    app.require_subcommand(1);

    CLI11_PARSE(app, argc, argv);
    cerr << "Command line options:" << full_cmd << endl;

    //run a sanity check on the input file provided
    if (!validate_gfa_input(inputpath)) {
        return 1;
    }

    if(*compact)run_compact(inputpath, outputpath, print_reverse, use_numeric, self_loop);
    else if(*decompose)run_decompose(inputpath, use_exact, minAlleles, iWalk);
}