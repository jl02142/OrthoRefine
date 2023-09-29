#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <filesystem> // std::filesystem::directory_iterator
#include <thread>
#include <future>
#include <math.h>     // abs overload
#include <cstring>    // std::strcpy

// The information from the feature table file is stored in the structure below
struct strct_feature_tables_info{           
    int size_strct;                         // size of structure created by new statement
    char chromosome = 'c';                  // can either be c (circular) or l (linear). Circular can overflow the ends when comparing gene's HOGs, linear cannot
    char domain = 'e';                      // b (bacteria), a (archaea) or e (eukaryote). bacteria and archaea have operons detected using the gene gap method if user turns option on
    std::vector<int> operon;                // Operon number as determined by gene gap method
    std::string name;                       // Feature table name
    std::vector<std::string> scaffold;      // Store scaffold number. Chromosome or complete level assembly have a single scaffold. Scaffolding and contigs assemblies have multiple scaffold numbers 
    std::vector<std::string> start_pos;     // Start position
    std::vector<std::string> end_pos;       // End position
    std::vector<std::string> prod_acc;      // product accession
    std::vector<std::string> gene_name;     // gene name
    std::vector<std::string> symbol;        // gene symbol
    std::vector<std::string> locus_tag;     // locus tag
    std::vector<char> direction;            // direction (+ or -)
    std::vector<int> HOG;                   // hierarchical orthogroups (HOGs) from OrthoFinder
};

namespace File_read{
    // read_ft reads and stores the feature table file information into the struct above
    std::vector<strct_feature_tables_info> * read_ft(std::string snf, std::string path, std::string OrthFnd_file, int diag, int ortho_finder_order);
    // count_HOG returns the number of lines in the OthoFinder file (N0.tsv). 
    int count_HOG(std::string path, std::string OrthFnd_file);
    // read_HOG reads the OrthoFinder file (N0.tsv) and calls assign_HOG and assign_HOG_match to match the prod_acc from the feature table to the 
    // OrthoFinder file (N0.tsv) and then loads the HOG information into hog_master.
    std::vector<int> *** read_HOG(std::string path, std::vector<strct_feature_tables_info> *feature_tables_info, int diag, std::string OrthFnd_file, int benchmark, int ortho_finder_order);
}

namespace File_handler{
    void is_number(std::string arg);
    bool file_exist(std::string file_name);     // verifies that a file exist for open and reading
    std::string load_ss(std::string file_name); // returns the contents of a file for processing, load string stream
}

// assign_HOG matches the prod_acc from the OrthoFinder file (N0.tsv) to the prod_acc of the feature table file. 
void assign_HOG(int temp_tab_count, int HOG_count, std::string temp_string, std::vector<strct_feature_tables_info> *feature_tables_info, std::vector<int> &i_storage, int benchmark, std::vector<std::string> name_storage, int ortho_finder_order);
// assign_HOG_match stores the location or oder of a gene/protein from the feature table into hog_master.
void assign_HOG_match(int numb_samp, int temp_tab_count, int HOG_count, int pairwise, std::vector<int> i_storage, std::vector<int>*** hog_master, std::vector<strct_feature_tables_info> *feature_tables_info, std::vector<std::string> name_storage, int ortho_finder_order);
// match_hog does the HOG comparison between the 2 gene windows
std::vector<std::vector<std::string>> match_hog(std::vector<int>*** HOG_master, std::vector<strct_feature_tables_info> *feature_tables_info, int pairwise, int HOG_line, int window_size, double synteny_ratio, int diag, int ft_gene_together, int print_all_orthofinder, int run_all_orthofinder, int benchmark, int run_single_HOG, int prod_acc_flag);
// while loop finds the bounds of the windows. Can change because of -999 (no data) or need to over or underflow in circular genomes
void while_loop(std::vector<int>*** HOG_master, int j, int p, int **zstore, int *z_neg_count, int *z_scaff_count, std::vector<strct_feature_tables_info> *feature_tables_info, int ftk, int m, int window_size, int diag, int hmnk, int pairwise, int combo[][2], int &z_max_operon_count);
// splits the locus_tag data structure in match_hog along tabs
void find_substr(bool &found_t, size_t &pos, std::vector<std::vector<std::string>> locus_tag, int xx, int yy, std::string &search);

int main(int argc, char* argv[]){

    std::string path = "";                  // path of files to be read; N0.tsv and feature table. 
    std::string refseq_file = "input.txt";  // user provided / created file that contains the REFSEQ accession per each sample. One accession per line in file
    std::string OrthFnd_file = "N0.tsv";    // Orthofinder file
    std::string outfile_name;               // outfile prefix name
    outfile_name = "outfile";
    int diag{0};                            // cout short diagnostic text if == 1 or long diagnostic text if == 2 or long long if == 3 or everything if == 4
    int window_size{8};                     // size to scan around a gene; this number does not include the gene itself in the count; window_size = 10 is 5 on each side
    double synteny_ratio{0.5};              // number of required matches within a window
    std::setprecision(2);
    int ft_gene_together{1};                // How should genes in a HOG be treated if they are next to eachother in the FT? 0 = print together if best; 1 = skip 
    int print_all_orthofinder{0};           // Controls printing all orthologs for benchmark file
    int benchmark{0};                       // Should the benchmark file be created?
    bool single_run{0};                     // changed to 1 when user provides window size or synteny ratio. If 0, run on programmed combinations of ws and sr
    int run_single_HOG{-1};                 // Run only on one HOG to save time
    int run_all_orthofinder{0};             // Run on all OrthoFinder results instead of those with at least one paralog
    int paralogs_print{0};                  // If paralogs should be considered for longest SOG printing. 0 = no, 1 = yes
    int prod_acc_flag{1};                   // If should print locus tag (1) or prod acc (0)
    int ortho_finder_order{0};              // If OrthoRefine output order should follow the input file generated by user (0) or follow the order from OrthoFinder's output

    for(int i = 1; i < argc; ++i){
        std::string arg = argv[i];
        if(arg == "-h" || arg == "--help"){
            std::cout << "Options are:" << '\n';
                std::cout << "--input is required and is the file created by the user that provides the REFSEQ accession per sample. One accession per line." << '\n';
                std::cout << "--OF_file is the OrthoFinder output file. Default is N0.tsv" << '\n';
                std::cout << "--window size controls how many total genes are considered when determining synteny for a single gene. if == 10, 5 genes on each side of the single gene being analyzed for synteny are used. Default is 8" << '\n';
                std::cout << "--synteny_ratio controls how many genes within a window must provide synteny support to classify the genes being compared as syntenous. Provide a deceimial point, i.e. 0.5 for half of the window. Default is 0.5." << '\n';
                std::cout << "--print_all controls if only HOGs with refinements are printed. 0 (only print those with changes and supported by synteny), 1 (print all including those without synteny support), or 2 (print all supported by syntenty, even if no changes)" << '\n';
                std::cout << "--run_all_orthofinder runs on all HOGs instead of only those with at least one paralog" << '\n';
                std::cout << "--outfile sets to the prefix to the outfile. Default name is outfile_ws_sr_pa_ra where ws is window size, sr is synteny ratio, pa is print all, and ra is run all" << '\n';
                std::cout << "--run_single_HOG will run only on the one HOG specified" << '\n';
                std::cout << "--path if input files are located at different location than current" << '\n';
                std::cout << "--diag will print extra information to diagnosis potential problems. Accepted values are 0 (none), 1 (short), or 2 (long) or 3 (long long) or 4 (everything). Default is 0" << '\n';
                std::cout << "--benchmark will print all ortholog pairs, confirmed by synteny, in a pairwise fashion (2 column output) for submission to the benchmarking web serivce" << '\n';
                std::cout << "--prod_acc controls if prod acc should be printed (0 default) or locus tag (1)" << '\n';
                //std::cout << "-together controls how the program handles when mutiple genes from a single sample are in a single HOG (and are thus being considered for synteny against the other members of the HOG from the other samples) and are next to each other in order in the feature table. Because there could be great difficultly in telling which gene is the 'better' ortholog (there may only be a differnce of one window match), we allow the program to consider all sequential genes from the HOG to be the best 'match' if any of the sequential genes are. Default is 1 and is recommended. 0 is do not allow the above to occur." << '\n';
            std::exit(EXIT_SUCCESS);
        }else if(arg == "--path"){
            path = argv[++i];
        }else if(arg == "--input"){
            refseq_file = argv[++i];
        }else if(arg == "--window_size"){
            File_handler::is_number(argv[++i]);
            window_size = std::stoi(argv[i]);
            single_run = 1;
        }else if(arg == "--OF_file"){
            OrthFnd_file.assign(argv[++i]);
        }else if(arg == "--synteny_ratio"){
            char *end_pointer;
            long int temp_numb = std::strtol(argv[++i], &end_pointer, 10);
            single_run = 1;
            if(end_pointer[0] == '.'){
                synteny_ratio = std::stof(argv[i]);
            }else{
                std::cout << "Error: synteny_ratio " << argv[i] << " is not in deceimial point." << '\n';
                std::exit(EXIT_FAILURE);
            }
        }else if(arg == "--paralogs_print"){
            paralogs_print = std::stoi(argv[++i]);
        }else if(arg == "--diag"){
            arg = argv[++i];
            if(arg == "0" || arg == "1" || arg == "2" || arg == "3" || arg == "4"){
                diag = std::stoi(argv[i]);
            }else{
                std::cout << "Error: -diag must be 0 (do not print diagnosis text) or 1 (print short text) or 2 (print long text) or 3(print long long text) or 4(print everything)." << '\n';
                std::exit(EXIT_FAILURE);
            }
        }else if(arg == "--together"){
            arg = argv[++i];
            if(arg == "0" || arg == "1"){
                ft_gene_together = std::stoi(argv[i]);
            }
        }else if(arg == "--ortho_finder_order"){
            ortho_finder_order = std::stoi(argv[++i]);
        }else if(arg == "--outfile"){
            outfile_name = argv[++i];
        }else if(arg == "--prod_acc"){
            prod_acc_flag = std::stoi(argv[++i]);
            if(prod_acc_flag != 0 && prod_acc_flag != 1){
                std::cout << "Error: prod_acc must be 0 (prod_acc) or 1 (locus_tag)" << '\n';
                std::exit(EXIT_FAILURE);
            }
        }else if(arg == "--print_all"){
            print_all_orthofinder = std::stoi(argv[++i]);
        }else if(arg == "--benchmark"){
            benchmark = std::stoi(argv[++i]);
        }else if(arg == "--run_single_HOG"){
            run_single_HOG = std::stoi(argv[++i]);
        }else if(arg == "--run_all_orthofinder"){
            run_all_orthofinder = std::stoi(argv[++i]);
        }else{
            std::cout << "Error: invalid argument provided." << '\t' << arg << '\n';
            std::exit(EXIT_FAILURE);
        }
    }            
    if(path == ""){ // If user does not provide a path
        path = std::filesystem::current_path();
    }

    // read the feature table files into a struct
    std::vector<strct_feature_tables_info> *feature_tables_info = File_read::read_ft(refseq_file, path, OrthFnd_file, diag, ortho_finder_order);
    // Number of lines in the "N0.tsv" encluding the header; number of total HOGS
    int HOG_line = File_read::count_HOG(path, OrthFnd_file);  
    // read OrthoFinder output (N0.tsv) and match the prod_acc from it to the prod_acc in the feature tables. 
    // Store the location / order from feature tables 
    std::vector<int> *** HOG_master = File_read::read_HOG(path, feature_tables_info, diag, OrthFnd_file, benchmark, ortho_finder_order);

    // pairwise = n(n - 1) / 2
    int pairwise = ((*feature_tables_info)[0].size_strct * ((*feature_tables_info)[0].size_strct - 1)) / 2;   

    // if == 0, run on predetermined combinations of window size and synteny ratio
    if(single_run == 0){
        int window_size_arr[10] = {2, 4, 6, 8, 10, 30, 40, 60, 80, 100}; // window size combinations
        double **synteny_ratio_arr = new double*[10];
        for(int wsa = 0; wsa < 10; ++wsa){  // wsa = window size array
            synteny_ratio_arr[wsa] = new double[4]{0};
            if(wsa < 4){
                // moved to below outside of wsa loop
            }else if(wsa < 8){
                double numb{0.5};
                for(int sra = 0; sra < 3; ++sra){  // synteny ratio array
                    synteny_ratio_arr[wsa][sra] = numb;
                    if(numb == 0.5){
                        numb = 0.3;
                    }else{
                        numb = 0.2;
                    }
                }
            }else{
                double numb{0.5};
                for(int sra = 0; sra < 4; ++sra){
                    synteny_ratio_arr[wsa][sra] = numb;
                    if(numb == 0.5){
                        numb = 0.3;
                    }else if(numb == 0.3){
                        numb = 0.2;
                    }else{
                        numb = 0.1;
                    }
                }
            }
        }
        synteny_ratio_arr[0][0] = 0.5;
        synteny_ratio_arr[1][0] = 0.25;
        synteny_ratio_arr[1][1] = 0.5;
        synteny_ratio_arr[2][0] = 0.2;
        synteny_ratio_arr[2][1] = 0.5;
        synteny_ratio_arr[3][0] = 0.2;
        synteny_ratio_arr[3][1] = 0.3;
        synteny_ratio_arr[3][2] = 0.5;

        /* // These are the combinations from above
        WS           |   SR
        2                0.5
        4                0.25 0.5
        6                0.20, 0.50
        8                0.20, 0.30, 0.50
        10, 20, 40, 60   0.20, 0.3, 0.5
        80, 100          0.10, 0.2, 0.3, 0.5
        */

        /* // print the combinations to verify
        for(int wsa = 0; wsa < 11; ++wsa){
            for(int sra = 0; sra < 4; ++sra){
                std::cout << wsa << '\t' << sra << '\t' << synteny_ratio_arr[wsa][sra] << '\n';
            }
            std::cout << '\n';
        }
        */

        std::vector<std::vector<std::vector<std::string>>> HOG_match;             // To store the result of the future
        std::vector<std::future<std::vector<std::vector<std::string>>>> threads;  // To store the thread calls in

        // Create thread calls and store them into vector called "threads". Loop loops through the possible window sizes and synteny ratios while calling function match_hog
        for(int wsa = 0; wsa < 10; ++wsa){ 
            window_size = window_size_arr[wsa];
            for(int sra = 0 ; sra < 4; ++sra){
                synteny_ratio = synteny_ratio_arr[wsa][sra];
                if(synteny_ratio == 0){
                    continue;
                }
                threads.push_back(std::async(std::launch::async, match_hog, HOG_master, feature_tables_info, pairwise, HOG_line, window_size, synteny_ratio, diag, ft_gene_together, print_all_orthofinder, run_all_orthofinder, benchmark, run_single_HOG, prod_acc_flag));
            }
        }

        // Loop through the vector called threads waiting for threads to finish before starting next thread. Store results into HOG_match. 
        for(auto &th : threads){
            HOG_match.push_back(th.get());
        }

        // final output statement
        std::vector<std::vector<std::vector<std::string>>> longest_chain(HOG_line, std::vector<std::vector<std::string>>(1));
        double ***record_ws_sr = new double**[HOG_line];      // array to store window size and synteny ratio values into
        double ***record_ws_sr_indx = new double**[HOG_line]; // array to store window size and synteny ratio index values into
        for(int i = 0; i < HOG_line; ++i){
            record_ws_sr[i] = new double*[100];
            record_ws_sr_indx[i] = new double*[100];
            for(int j = 0; j < 100; ++j){
                record_ws_sr[i][j] = new double[2];
                record_ws_sr_indx[i][j] = new double[2];
                for(int k = 0; k < 2; ++k){
                    record_ws_sr[i][j][k] = 0;
                    record_ws_sr_indx[i][j][k] = 0;
                }
            }
        }

        int ws = 0;      // window size
        int ws_indx{0};  // window size index
        double sr = 0;   // synteny ratio
        int sr_indx = 0; // synteny_ratio index
        double sog_per_combo[HOG_match.size()];
        for(int i = 0; i < HOG_match.size(); ++i){
            sog_per_combo[i] = 0;
        }
        double jan_avg[HOG_match.size()]{0};
        int longest_chain_count[HOG_line]{0};
        for(int h = 0; h < HOG_match.size(); ++h){ // loop through the results
            ws = window_size_arr[ws_indx];
            sr = synteny_ratio_arr[ws_indx][sr_indx];
            while(sr == 0){
                if(sr_indx == 3){
                    ++ws_indx;
                    sr_indx = 0;
                }else{
                    ++sr_indx;
                }
                ws = window_size_arr[ws_indx];
                sr = synteny_ratio_arr[ws_indx][sr_indx];
            }

            // https://stackoverflow.com/questions/57882748/remove-trailing-zero-in-c
            std::cout << ws << '\t' << sr << '\t' << ws_indx << '\t' << sr_indx << '\t' << "--------------------------------------------------------------------" << '\n';
            std::stringstream ss;
            ss <<  std::to_string(sr);
            std::string temp = ss.str();
            temp = temp.substr(0, temp.find_last_not_of('0')+1); // strip extra 0's
            std::string outfile_name_2 = outfile_name + '_' + std::to_string(ws) + '_' + temp + '_' + std::to_string(print_all_orthofinder) + '_' + std::to_string(run_all_orthofinder);
            std::ofstream outfile(outfile_name_2);

            int final_refine{0};      // total number of refinements, HOGs may be counted twice if a HOG can be split mutiple times.
            int total_hog_refine{0};  // total number of HOGs refined. HOGs will only be counted once regardless of how many splits. 
            bool flag_hog_counted{0}; // keep track is a hog is already counted for total_hog_refine

            int hog_count{0};
            for(int i = 0; i < HOG_line; ++i){                          // loop through each HOG
            //if(i != 2){
            //    continue;
            //}
                int within_longest_chain_count{0};
                int new_long_chain{0};                                  // 0 no new long chain, 1 is new long chain, 2 is tie
                if(HOG_match[h][i][0] != ""){
                    ++hog_count;
                }
                std::vector<int> remember_j;

                outfile << "HOG" << '\t' << "SOG" << '\t' << "Gene_name" << '\t'; // outfile header line
                for(int i = 0; i < (*feature_tables_info)[0].size_strct; ++i){
                    outfile << (*feature_tables_info)[i].name;
                    if(i < (*feature_tables_info)[0].size_strct - 1){
                        outfile << '\t';
                    }
                }
                outfile << '\n';
                int sog_count{0};
                for(int j = 0; j < HOG_match[h][i].size(); ++j){        // loop through each SOG
                    int chain_tab_count{0};

                    for(int k = 0; k < HOG_match[h][i][j].size(); ++k){ // loop through each member of the SOG
                        if(HOG_match[h][i][j][k] == '\t'){              // count number of tabs (genomes) to find longest SOG
                            ++chain_tab_count;
                        }else if(HOG_match[h][i][j][k] == ','){         // if contains a paralog (contains a ","), don't allow to be considered as longest
                            if(paralogs_print == 0)goto skip_chain;
                        }
                    }
                    if(chain_tab_count > longest_chain_count[i]){ // new longest SOG
                        longest_chain_count[i] = chain_tab_count; // store value (number of genomes) of longest SOG 
                        new_long_chain = 1;
                        remember_j.clear();
                        remember_j.push_back(j);
                        if(diag > 2)std::cout << i << '\t' << "YES NEW LONGEST CHAIN" << '\t' << chain_tab_count << '\t' << longest_chain_count[i] << '\t' << "ws" << '\t' << ws << '\t' << "sr" << '\t' << sr << '\n';
                    }else if(chain_tab_count > 0 && chain_tab_count == longest_chain_count[i]){ // tie in longest SOG
                        if(new_long_chain == 0)new_long_chain = 2;
                        remember_j.push_back(j);
                        if(diag > 2)std::cout << i << '\t' << "SAME MATCH CHAIN" << '\t' << chain_tab_count << '\t' << longest_chain_count[i] << '\t' << "ws" << '\t' << ws << '\t' << "sr" << '\t' << sr << '\n';
                    }else{ // No new or tie in longest SOG
                        if(diag > 2)std::cout << i << '\t' << "NO NEW LONGEST CHAIN" << '\t' << chain_tab_count << '\t' << longest_chain_count[i] << '\t' << "ws" << '\t' << ws << '\t' << "sr" << '\t' << sr <<  '\n';
                    }
                    if(chain_tab_count > within_longest_chain_count){
                        within_longest_chain_count = chain_tab_count;
                    }
                    skip_chain:;
                    if(HOG_match[h][i][j] != ""){
                        if(flag_hog_counted == 0){
                            ++total_hog_refine;
                            flag_hog_counted = 1;
                        }
                        ++final_refine;
                        outfile << "N0.HOG" << std::setfill('0') << std::setw(7) << i << '\t' << i << '.'<<  sog_count << '\t' <<
                        (*feature_tables_info)[0].gene_name[HOG_master[0][i][0][0]] << '\t' << HOG_match[h][i][j] << '\n';
                        ++sog_count;
                    }
                
                }
                
                if(new_long_chain == 1){  // new longest so replace
                    longest_chain[i][0].erase(longest_chain[i][0].begin(), longest_chain[i][0].end());
                    for(int j = 0; j < remember_j.size(); ++j){
                        longest_chain[i][0].push_back(HOG_match[h][i][remember_j[j]]);
                        record_ws_sr[i][j][0] = ws;
                        record_ws_sr[i][j][1] = sr;
                        record_ws_sr_indx[i][j][0] = ws_indx;
                        record_ws_sr_indx[i][j][1] = sr_indx;
                    }
                }else if(new_long_chain == 2){ // tie
                    for(int j = 0; j < remember_j.size(); ++j){
                        longest_chain[i][0].push_back(HOG_match[h][i][remember_j[j]]);
                        int k = longest_chain[i][0].size() - 1; // can prob just use m
                        record_ws_sr[i][k][0] = ws;
                        record_ws_sr[i][k][1] = sr;
                        record_ws_sr_indx[i][k][0] = ws_indx;
                        record_ws_sr_indx[i][k][1] = sr_indx;
                    }
                }
                flag_hog_counted = 0;
                sog_per_combo[h] += within_longest_chain_count; // for calculating average below
            }
            outfile << "Number of HOGs refined:" << '\t' << total_hog_refine << '\t' << "for a total refinement of" << '\t' << final_refine << '\n';
            outfile.close();
            if(benchmark == 1){ // delete output file as a blank file was created. 
                char *cstring = new char[outfile_name_2.length() +1];
                std::strcpy(cstring, outfile_name_2.c_str());
                std::remove(cstring);
            }
            if(sr_indx == 3){
                ++ws_indx;
                sr_indx = 0;
            }else{
                ++sr_indx;
            }
            //std::cout << "AVG" << '\t' << double(sog_per_combo[h]) / hog_count << '\t' << sog_per_combo[h] << '\t' << hog_count<< '\n';
            jan_avg[h] = double(sog_per_combo[h]) / hog_count;
        }
        
        std::string outfile_name_3 = "longest_SOG.tsv";
        std::ofstream outfile3(outfile_name_3);
        for(int i = 0; i < HOG_line; ++i){
            int j = 0;
                for(int k = 0; k < longest_chain[i][j].size(); ++k){
                    if(longest_chain[i][j][k].size() < 1){
                        continue;
                    }
                    int count{0};
                    for(int m = 0; m < longest_chain[i][j][k].size(); ++m){
                        if(longest_chain[i][j][k][m] == '\t'){
                            ++count;
                        }
                    }
                    outfile3 << "HOG" << '\t' << i << '\t' << j << '\t' << "ws" << '\t' << record_ws_sr[i][k][0] << '\t' << "sr" << '\t' << record_ws_sr[i][k][1] << '\t' << count << '\t' << longest_chain[i][j][k] <<  '\n';
                }
        }
        outfile3.close();
        if(benchmark == 1){ // delete longest sog file as a blank file was created. 
            char *cstring = new char[outfile_name_3.length() +1];
            std::strcpy(cstring, outfile_name_3.c_str());
            std::remove(cstring);
        }
                
        // print the average. default is longest sog with no paralogs, can be changed to longest sog with paralogs with paralog option
        int need_count{0};
        for(int i = 0; i < 10; ++i){
            int ws = window_size_arr[i];
            for(int j = 0; j < 4; ++j){
                double sr = synteny_ratio_arr[i][j];
                if(sr == 0){
                    continue;
                }
                std::cout << "AVG" << '\t' << ws << '\t' << sr << '\t' << jan_avg[need_count] << '\n';
                ++need_count;
            }
        }
    }else{  // for running on user provided window_size and synteny_ratio
        // match_hog performs the synteny analysis in the sliding window 
        std::vector<std::vector<std::string>> HOG_match = match_hog(HOG_master, feature_tables_info, pairwise, HOG_line, window_size, synteny_ratio, diag, ft_gene_together, print_all_orthofinder, run_all_orthofinder, benchmark, run_single_HOG, prod_acc_flag);

        // final output statement
        std::cout << window_size << '\t' << synteny_ratio << '\t' << "--------------------------------------------------------------------" << '\n';
        std::stringstream ss;
        ss <<  std::to_string(synteny_ratio);
        std::string temp = ss.str();
        temp = temp.substr(0, temp.find_last_not_of('0')+1);  // remove extra 0's
        std::string outfile_name_2 = outfile_name + '_' + std::to_string(window_size) + '_' + temp + '_' + std::to_string(print_all_orthofinder) + '_' + std::to_string(run_all_orthofinder);
        std::ofstream outfile(outfile_name_2);

        int final_refine{0};      // total number of refinements, HOGs may be counted twice if a HOG can be split mutiple times.
        int total_hog_refine{0};  // total number of HOGs refined. HOGs will only be counted once regardless of how many splits. 
        bool flag_hog_counted{0}; // keep track is a hog is already counted for total_hog_refine

        outfile << "HOG" << '\t' << "SOG" << '\t' << "Gene_name" << '\t'; // outfile header line
        for(int i = 0; i < (*feature_tables_info)[0].size_strct; ++i){
            outfile << (*feature_tables_info)[i].name;
            if(i < (*feature_tables_info)[0].size_strct - 1){
                outfile << '\t';
            }
        }
        outfile << '\n';

        for(int i = 0; i < HOG_line; ++i){
            int sog_count{0};
            for(int j = 0; j < HOG_match[i].size(); ++j){
                if(HOG_match[i][j] != ""){
                    if(flag_hog_counted == 0){
                        ++total_hog_refine;
                        flag_hog_counted = 1;
                    }
                    ++final_refine;
                    outfile << "N0.HOG" << std::setfill('0') << std::setw(7) << i << '\t' << i << '.'<<  sog_count << '\t' << (*feature_tables_info)[0].gene_name[HOG_master[0][i][0][0]] << '\t' << HOG_match[i][j] << '\n';
                    ++sog_count;
                }
            }
            flag_hog_counted = 0;
        }
        outfile << "Number of HOGs refined:" << '\t' << total_hog_refine << '\t' << "for a total refinement of" << '\t' << final_refine << '\n';
        outfile.close();
        if(benchmark == 1){ // delete output file as blank file was created
            char *cstring = new char[outfile_name_2.length() +1];
            std::strcpy(cstring, outfile_name_2.c_str());
            std::remove(cstring);
        }
    }
    std::cout << "Done" << '\n';
}

namespace File_handler{
    void is_number(std::string arg){
        for(int i = 0; i < arg.length(); ++i){
            if(isdigit(arg[i]) == false){
                std::cout << "Error : " << arg << " is not an integer." << std::endl;
                std::exit(EXIT_FAILURE);
            }
        }
    }
}

namespace File_handler{
    bool file_exist(std::string file_name){
        std::ifstream file(file_name);
        if(file){
            return 1;
        }else{
            std::cout << "Error: file " << file_name << " not found" << "\n";
            std::exit(EXIT_FAILURE);
        }
    }
}

namespace File_handler{
    std::string load_ss(std::string file_name){
        File_handler::file_exist(file_name);
        std::string temp_string;
        std::ifstream infile(file_name);
        infile.seekg(0, std::ios::end);     // Go to end of file
        temp_string.resize(infile.tellg()); // resize temp_string to size of infile
        infile.seekg(0, std::ios::beg);     // Go to start of file
        infile.read(&temp_string[0], temp_string.size());  // Read file into temp_string
        infile.close();
        return temp_string;
    }
}


namespace File_read{
    // read the info from each feature table
    std::vector<strct_feature_tables_info> * read_ft(std::string snf, std::string path, std::string OrthFnd_file, int diag, int ortho_finder_order){ // snf = sample name file
        
        int path_length = path.length() + 1;
        std::stringstream ss;
        char c;
        std::string temp_string = "";
        int temp_tab_count{0};
        bool second_col{0};  // second column from user created input file, checks if user provided the second column
        bool third_col{0};

        ss << (File_handler::load_ss(snf));

        while(ss >> std::noskipws >> c){  // reads user generated input file with GCF accession to count how many lines and if second and third column are present
            if(c == '\n'){
                ++temp_tab_count; 
            }else if(c == '\t' && second_col == 0){
                second_col = 1;
            }else if(c == '\t' && second_col == 1){
                third_col = 1;
            }
        }
        if(second_col == 1 && third_col == 0){
            std::cout << "Error: must specify both chromosome shape: linear (l) or circular (c) in second column AND domain: bacteria (b), Archaea (a), or Eukarya (e) in third column" << '\n';
            std::exit(EXIT_FAILURE);
        }
        std::string GCF_prefix[temp_tab_count]; // Array to store the GCF accession 
        int GCF_length = temp_tab_count;        // Equal to the number of input samples listed in input file by user
        temp_tab_count = 0;
        ss.clear();
        ss.seekg(0);
        while(ss >> std::noskipws >> c){        // Store the GCF_prefix provided by the user
            if(c == '\n' || c == '\t'){
                if(temp_string.size() > 1){
                    GCF_prefix[temp_tab_count] = temp_string;
                    ++temp_tab_count; 
                }else if(temp_string != "c" && temp_string != "C" && temp_string != "l" && temp_string != "L" && temp_string != "a" && temp_string != "b" && temp_string != "e" && temp_string != "A" && temp_string != "B" && temp_string != "E"){
                    std::cout << temp_string << '\n';
                    std::cout << "ERROR: second column contains non 'c' 'C' 'l' 'L'" << '\n';
                    std::exit(EXIT_FAILURE);
                }
                temp_string = "";
            }else if(c == ' '){
                std::cout << "Error: Space detected when tab is needed, probably between first and second column on line" << '\t' << temp_tab_count << '\n';
                std::exit(EXIT_FAILURE);
            }else{
                temp_string += c;
            }            
        }

        ss.clear();
        ss.str("");
        bool OF_file_found = 0; // OrthoFinder file
        for(const auto &file : std::filesystem::directory_iterator(path)){ // loops through each file name in current dir to match to provided name by user at run time for OrthoFinder output
            std::string file_path = file.path();
            std::string compare = file_path.substr(path_length);  // Remove full path; keep just file name
            std::size_t found = compare.find(OrthFnd_file);
            if(found != std::string::npos && OF_file_found == 1){ // can prob remove OF_file_found flag as we stopped using recursive::directory_iterator and just use directory_iterator
                std::cout << "ERROR: multiple matching OrthoFinder files named " << OrthFnd_file << " found in this directory or all sub-directories. Change duplicate file names." << '\n';
                std::exit(EXIT_FAILURE);
            }
            if(found != std::string::npos){
                ss << File_handler::load_ss(OrthFnd_file);
                OF_file_found = 1;
            }            
        }

        temp_string = "";
        temp_tab_count = 0;
        while(ss >> std::noskipws >> c){  // Count number of tabs present in the file's first line
            if(c == '\t'){
                ++temp_tab_count;
            }else if(c == '\n'){
                break;
            }
        }
        std::string ordered_GCF_input[temp_tab_count];
        temp_tab_count = 0;
        ss.seekg(0); // Go back to start of file

        while(ss >> std::noskipws >> c){  // Store the GCF prefix in the same order as the OrthFinder output
            if(c != '\t' && c != '\n'){
                temp_string += c;
            }else if(ortho_finder_order == 0){
                if(c == '\t' || c == '\n'){
                    if(temp_tab_count >= 3){
                        for(int i = 0; i < GCF_length; ++i){
                            std::string GCF_name_temp = GCF_prefix[i];
                            if(temp_string.find(GCF_name_temp) == 0){
                                ordered_GCF_input[i] = temp_string;
                            }
                        }
                    }
                    temp_string = "";
                    ++temp_tab_count;
                    if(c == '\n'){
                        break;
                    }
                }
            }else if(c == '\t'){
                if(temp_tab_count >= 3){
                    for(int i = 0; i < GCF_length; ++i){
                        std::string GCF_name_temp = GCF_prefix[i];
                        if(temp_string.find(GCF_name_temp) == 0){
                            ordered_GCF_input[temp_tab_count - 3] = temp_string;
                        }
                    }

                }
                temp_string = "";
                ++temp_tab_count;
            }else if(c == '\n'){
                for(int i = 0; i < GCF_length; ++i){
                    std::string GCF_name_temp = GCF_prefix[i];
                    if(temp_string.find(GCF_name_temp) == 0){
                        ordered_GCF_input[temp_tab_count - 3] = temp_string;
                    }
                }                
                break;
            }
        }

        for(int i = 0; i < GCF_length; ++i){
            std::cout << GCF_prefix[i] << '\n';
        }
        // Uncomment below loop to verify order GCF prefix
        for(std::string i : ordered_GCF_input){
            std::cout << i << '\t';
        }
        std::cout << '\n';
        
        std::vector<std::string> ft_file_names(GCF_length);
        ss.clear(); 
        ss.str(""); // clear previous file "N0.tsv" from ss
        ss << (File_handler::load_ss(snf));
        temp_string = "";

        std::vector<strct_feature_tables_info> *feature_tables_info = new std::vector<strct_feature_tables_info>(GCF_length);

        int remember_i{-1};
        temp_tab_count = 0;
        // This loop is not the best when snf is already in memory but I'm patching in reading snf (above code) afterwards
        while(ss >> std::noskipws >> c){
            if(c == '\n' || c == '\t'){
                bool file_match{0};
                if(c == '\t' && second_col == 1){ // if second column is present in input file, skip over checking if the feature table exist based on the single character in the second col. always would fail
                    (*feature_tables_info)[remember_i].chromosome = temp_string[0];  // store the char from the second col
                    file_match = 1;
                }else if(c == '\n' && third_col == 1){
                    (*feature_tables_info)[remember_i].domain = temp_string[0]; // store the char from the third col
                    file_match = 1;
                }else{

                }
                if(temp_tab_count == 0){
                    for(const auto &file : std::filesystem::directory_iterator(path)){
                        std::string file_path = file.path();
                        std::string compare = file_path.substr(path_length);  // Remove full path; keep just file name
                        std::size_t found = compare.find(temp_string);  // Compare file name in dir to file name provided by user
                        if(found != std::string::npos){
                            if(compare.find("ft_locus_tag_replaced_with_uniprot") != std::string::npos){  // used when benchmarking
                                found = compare.find("ft_locus_tag_replaced_with_uniprot");
                                std::cout << "Warning: Using modified feature table file."  << '\t' <<  compare << '\n';
                            }else{
                                found = compare.find("feature_table.txt");  // Check if file name in dir also has "feature_table" in name
                            }
                            if(found != std::string::npos){
                                for(int i = 0; i < GCF_length; ++i){
                                    
                                    std::string GCF_name_temp = ordered_GCF_input[i]; //GCF_name_temp
                                    for(int j = 0; j < GCF_length; ++j){
                                        std::string ordered_GCF_name_temp = GCF_prefix[j]; //ordered_GCF_name_temp
                                        std::size_t pos = GCF_name_temp.find(ordered_GCF_name_temp);
                                        if(pos != std::string::npos){
                                            if(compare.find(ordered_GCF_name_temp) == 0){
                                                std::cout << compare << '\n';
                                                ft_file_names[i] = compare;
                                                remember_i = i;
                                                file_match = 1;
                                                goto break_auto;
                                            }else{

                                            }
                                        }
                                    }                            
                                }
                            }
                        }
                    }
                }
                break_auto:;

                if(file_match == 0){
                    std::cout << "Error feature table file missing" << '\n';
                    std::exit(EXIT_FAILURE);
                }
                temp_string = "";
                if(c == '\t'){
                    ++temp_tab_count;
                }else{ // is \n
                    temp_tab_count = 0;
                }
            }else{
                temp_string += c;
            }
        }

        /* // Uncomment below loop to print feature table file names stored in memory
        std::cout << "PRINT FT FILE NAMES" << '\n';
        for(std::string i : ft_file_names){
            std::cout << i << '\n';
        }
        std::cout << "END PRINT FT FILE NAMES" << '\n';
        */

        int ft_line_count{0};
        temp_tab_count = 0;
        bool h_flag{0}; // header line flag
        bool c_flag{0}; // line starts with "CDS"
        for(int i = 0; i < ft_file_names.size(); ++i){

            ft_line_count = 0;
            (*feature_tables_info)[i].name = ft_file_names[i];
            ss.clear();
            ss << File_handler::load_ss(ft_file_names[i]);
            while(ss >> std::noskipws >> c){  // Count how many lines in feature table are "CDS" but also not "psuedogene".
                if(h_flag == 0){
                    if(c == '\n'){
                        h_flag = 1; // Header line
                    }
                }else if(h_flag == 1){
                    if(c == '\t'){
                        ++temp_tab_count;
                    }else if(temp_tab_count == 0 && c == 'C'){ // Only count lines start with "CDS"
                        c_flag = 1;
                        ++ft_line_count;
                    }else if(c_flag == 1 && temp_tab_count == 1 && c == 'u'){
                        --ft_line_count; // Don't want psuedo genes, match to 'u' in 'without'
                    }else if(c == '\n'){
                        c_flag = 0;
                        temp_tab_count = 0;
                    }
                }
            }

            (*feature_tables_info)[i].size_strct = ft_file_names.size();
            (*feature_tables_info)[i].scaffold.resize(ft_line_count);
            (*feature_tables_info)[i].start_pos.resize(ft_line_count);
            (*feature_tables_info)[i].end_pos.resize(ft_line_count);
            (*feature_tables_info)[i].prod_acc.resize(ft_line_count);
            (*feature_tables_info)[i].gene_name.resize(ft_line_count);
            (*feature_tables_info)[i].symbol.resize(ft_line_count);
            (*feature_tables_info)[i].locus_tag.resize(ft_line_count);
            (*feature_tables_info)[i].HOG.resize(ft_line_count);
            (*feature_tables_info)[i].operon.resize(ft_line_count);
            (*feature_tables_info)[i].direction.resize(ft_line_count);
            for(int j = 0; j < ft_line_count; ++j){
                (*feature_tables_info)[i].scaffold[j] = "";
                (*feature_tables_info)[i].start_pos[j] = "";
                (*feature_tables_info)[i].end_pos[j] = "";
                (*feature_tables_info)[i].prod_acc[j] = "";
                (*feature_tables_info)[i].gene_name[j] = "";
                (*feature_tables_info)[i].symbol[j] = "";
                (*feature_tables_info)[i].locus_tag[j] = "";
                (*feature_tables_info)[i].direction[j] = '\0';
                (*feature_tables_info)[i].HOG[j] = -999;  // Store a number that we will never observe in real data so we know when something funky is up
                (*feature_tables_info)[i].operon[j] = -1;
            }
            if(diag > 2){
                std::cout << "feature_tables_name" << '\t' << (*feature_tables_info)[i].name << '\t' << "AND SIZE" << '\t' << ft_line_count << '\n';
            }

            ft_line_count = 0;
            temp_tab_count = 0;
            h_flag = 0;
            c_flag = 0;
            ss.clear();
            ss << File_handler::load_ss(ft_file_names[i]);
            int sspos;  // To store location in file where line starts. stringstream position     
            while(ss >> std::noskipws >> c){
                
                if(h_flag == 0){
                    if(c == '\n'){
                        h_flag = 1;
                    }
                }else if(h_flag == 1){
                    if(c == '\t'){
                        ++temp_tab_count;
                    }else if(temp_tab_count == 0 && c == 'C'){  // Only lines that start with "CDS"
                        sspos = ss.tellg();
                        c_flag = 1;
                    }else if(c_flag == 1 && temp_tab_count == 1 && c == 'u'){  // Don't want psuedogenes, look for "u" in "without"
                        c_flag = 0;
                    }

                    if(c_flag == 1 && temp_tab_count == 2){ // Need any tab # after 1 to allow skipping of psuedogenes
                        temp_tab_count = 0;
                        ss.seekg(sspos - 1);  // Return to start of line
                        bool tab_6_flag{0};
                        std::string scaffold_string = "";
                        while(ss >> std::noskipws >> c){
                            if(c == '\t'){
                                ++temp_tab_count;
                                continue;
                            }
                            if(temp_tab_count == 6){ // record scaffold name
                                scaffold_string += c;
                                tab_6_flag = 1;
                            }else if(temp_tab_count == 7){  // Start pos
                                (*feature_tables_info)[i].start_pos[ft_line_count] += c;
                                if(tab_6_flag == 1){
                                    (*feature_tables_info)[i].scaffold[ft_line_count] = scaffold_string;
                                    scaffold_string = "";
                                    tab_6_flag = 0;
                                }
                            }else if(temp_tab_count == 8){  // End pos
                                (*feature_tables_info)[i].end_pos[ft_line_count] += c;
                            }else if(temp_tab_count == 9){
                                (*feature_tables_info)[i].direction[ft_line_count] = c;
                            }else if(temp_tab_count == 10){  // Product accession
                                (*feature_tables_info)[i].prod_acc[ft_line_count] += c;
                            }else if(temp_tab_count == 13){  // Name 
                                (*feature_tables_info)[i].gene_name[ft_line_count] += c;
                            }else if(temp_tab_count == 14){  // Symbol
                                (*feature_tables_info)[i].symbol[ft_line_count] += c;
                            }else if(temp_tab_count == 16){ // Locus tag
                                (*feature_tables_info)[i].locus_tag[ft_line_count] += c;
                            }
                            if(c == '\n'){
                                ++ft_line_count;
                                c_flag = 0;
                                temp_tab_count = 0;
                                break;
                            }
                        }
                    }
                }
                if(c == '\n'){  // Reset after non "CDS" line
                    temp_tab_count = 0;
                }
            }
            
            //std::cout << "ft_line_count" << '\t' << ft_line_count << '\n';
            //std::cout << "name" << '\t' << (*feature_tables_info)[i].name << '\n';
            int operon_count{0};
            if((*feature_tables_info)[i].domain == 'b' || (*feature_tables_info)[i].domain == 'a'){
                for(int j = 0; j < ft_line_count - 1; ++j){
                    if((*feature_tables_info)[i].scaffold[j] == (*feature_tables_info)[i].scaffold[j + 1]){
                        // Gene gap method from Yan and Moult 2006. 
                        if(abs(std::stoi((*feature_tables_info)[i].start_pos[j + 1]) - std::stoi((*feature_tables_info)[i].end_pos [j])) < 25){
                            if((*feature_tables_info)[i].direction[j] == (*feature_tables_info)[i].direction[j + 1]){
                                //std::cout << "POS" << '\t' << (*feature_tables_info)[i].start_pos[j + 1] << '\t' << (*feature_tables_info)[i].end_pos[j] << '\t' << (*feature_tables_info)[i].  direction[j + 1] << '\t' << (*feature_tables_info)[i].direction[j] << '\n';
                                (*feature_tables_info)[i].operon[j] = operon_count;
                                (*feature_tables_info)[i].operon[j + 1] = operon_count;
                            }else{
                                //std::cout << "NOT direction" << '\t' << (*feature_tables_info)[i].start_pos[j + 1] << '\t' << (*feature_tables_info)[i].end_pos [j] << '\t' <<  (*feature_tables_info)[i].direction[j + 1] << '\t' << (*feature_tables_info)[i].direction[j] << '\n';
                            }
                        }else{
                            //std::cout << "NOT POS" << '\t' << (*feature_tables_info)[i].start_pos[j + 1] << '\t' << (*feature_tables_info)[i].end_pos [j] << '\n';
                            ++operon_count; // this counter can count mutiple single gene operons, still works but does not report true number of operons. which we don't need
                        }
                    }else{
                        //std::cout << "NOT SCAFFOLD" << '\t' << j << '\t' << (*feature_tables_info)[i].scaffold[j] << '\t' << (*feature_tables_info)[i].scaffold[j + 1] << '\n';
                    }
                }
                for(int j = 0; j < ft_line_count; ++j){
                    //std::cout << "operon" << '\t' << j << '\t' << (*feature_tables_info)[i].operon[j] << '\t' << (*feature_tables_info)[i].start_pos[j] << '\n';
                }
            }
        }
        return feature_tables_info;
    }

}

namespace File_read{  
    // code repeats partly from read_HOG but this was easier way to return file line count
    int count_HOG(std::string path, std::string OrthFnd_file){
        File_handler::file_exist(OrthFnd_file); 
        std::stringstream ss;
        int path_length = path.length();
        for(const auto &file : std::filesystem::directory_iterator(path)){ // was recursive_directory_iterator
            std::string file_path = file.path();
            std::string compare = file_path.substr(path_length);
            std::size_t found = compare.find(OrthFnd_file);
            if(found != std::string::npos){
                ss << File_handler::load_ss(file_path);
                path = file.path();
                break;
            }  
        }
        char c;
        bool h_flag{0};  // Flag to skip reading in to memory header line
        int hog_line_count{0};
        while(ss >> std::noskipws >> c){ // count number of lines in OrthoFinder output file
            if(h_flag == 0){
                if(c == '\n'){
                    h_flag = 1;
                }
            }else if(h_flag == 1){
                if(c == '\n'){
                    ++hog_line_count;
                }
            }
        }
        return hog_line_count;
    }
}

namespace File_read{

    std::vector<int> *** read_HOG(std::string path, std::vector<strct_feature_tables_info> *feature_tables_info, int diag, std::string OrthFnd_file, int benchmark, int ortho_finder_order){

        std::stringstream ss;
        int path_length = path.length();
        for(const auto &file : std::filesystem::directory_iterator(path)){
            std::string file_path = file.path();
            std::string compare = file_path.substr(path_length);
            std::size_t found = compare.find(OrthFnd_file);
            if(found != std::string::npos){
                ss << File_handler::load_ss(file_path);
                path = file.path();  
                break;
            }  
        }

        char c;
        bool h_flag{0};  // Flag to skip reading in to memory header line
        int hog_line_count{0};
        while(ss >> std::noskipws >> c){  // count number of lines in OrthoFinder output file
            if(h_flag == 0){
                if(c == '\n'){
                    h_flag = 1;
                }else if(c == '\r'){
                    std::cout << "ERROR: N0.tsv file appears to be in windows format with \\r\\n as newline character. Try:" << std::endl;
                    std::cout << "dos2unix N0.tsv" << '\n';
                    std::cout << "And then rerun OrthoRefine" << '\n';
                    std::exit(EXIT_FAILURE);
                }
            }else if(h_flag == 1){
                if(c == '\n'){
                    ++hog_line_count;
                }
            }
        }
        /*
         // Uncomment to see number of lines in OrthoFinder output file
        std::cout << "HOG LINE COUNT" << '\t' << hog_line_count << '\n';
        */

        int p = (*feature_tables_info)[0].size_strct;  // number of samples
        int pairwise = (p * (p - 1)) / 2;

        std::vector<int> ***hog_master = new std::vector<int>**[pairwise];
        for(int i = 0; i < pairwise; ++i){
            hog_master[i] = new std::vector<int>*[hog_line_count];
            for(int j = 0; j < hog_line_count; ++j){
                hog_master[i][j] = new std::vector<int>[2];
                for(int k = 0; k < 2; ++k){
                    hog_master[i][j][k].resize(1);
                    hog_master[i][j][k][0] = -9999;
                }
            }
        }

        // hog_master details
        // First level or "i" is the pairwise number
        // Second level or "j" is the HOG number from OrthoFinder's output
        // Third level or "k" is which of the two samples from the pairwise comparison
        // Fourth level or "0" is for mutiple matches from OrthoFinder's output per sample per HOG

        ss.clear();
        ss << File_handler::load_ss(path);
        h_flag = 0;            // Header line flag
        bool double_tab_f{0};  // Double tab flag to detect a double tab (no match for sample to HOG) in the NO.tsv file
        int h{-1};
        int temp_tab_count{0};
        int file_tab_count{1}; // Total number of tabs on one line
        int HOG_count{0};
        int numb_samp = p;     // number of samples provided by user
        std::string temp_string = "";
        std::vector<int> i_storage;  // Vector to store what location / order a gene is in the feature table
        std::vector<std::string> name_storage; // Vector to store names of GCF files from N0.tsv 
        
        while(ss >> std::noskipws >> c){  // Match prod_acc from OrthoFinder's output and feature tables. Store gene order from feature tables in hog_master.

            if(h_flag == 0){
                if(c == '\t'){
                    if(file_tab_count > 3){
                        name_storage.push_back(temp_string);
                    }
                    ++file_tab_count;
                    temp_string = "";  
                }else if(c == '\n'){
                    h_flag = 1;  // Header flag to skip the first line
                    name_storage.push_back(temp_string);
                    temp_string = "";
                }else{
                    temp_string += c;
                }
            }else if(h_flag == 1){
                if(c != '\t'){
                    double_tab_f = 0;  // two tabs can happen when a sample does not provide a gene to a particular HOG as tabs are used to delimit between samples by OrthoFinder
                }
                if(c == '\t'){
                    if(double_tab_f == 1){ // check for missing data; sample did not contrinbute a gene to HOG
                        ++temp_tab_count;
                        continue;
                    }
                    if(temp_tab_count >= 3){
                        assign_HOG(temp_tab_count, HOG_count, temp_string, feature_tables_info, i_storage, benchmark, name_storage, ortho_finder_order);  // matches prod_acc from OrthoFinder to prod_acc in feature table
                        assign_HOG_match(numb_samp, temp_tab_count, HOG_count, pairwise, i_storage, hog_master, feature_tables_info, name_storage, ortho_finder_order);  // loads into hog_master the location / order of the prod_acc from the feature table
                        i_storage.clear(); 
                    }
                    ++temp_tab_count;
                    temp_string = "";
                    double_tab_f = 1;
                }else if(c == ','){  // OrthoFinder uses comma to delimit when multiple genes from a sample are members of a single HOG. Probably paralog present.
                    assign_HOG(temp_tab_count, HOG_count, temp_string, feature_tables_info, i_storage, benchmark, name_storage, ortho_finder_order);                    
                    temp_string = "";
                }else if(c == '\n'){
                    if(temp_string == ""){ // if final sample did not have a member for that hog
                        temp_tab_count = 0;
                        ++HOG_count;
                        i_storage.clear();
                        continue;
                    }
                    assign_HOG(temp_tab_count, HOG_count, temp_string, feature_tables_info, i_storage, benchmark, name_storage, ortho_finder_order);
                    assign_HOG_match(numb_samp, temp_tab_count, HOG_count, pairwise, i_storage, hog_master, feature_tables_info, name_storage, ortho_finder_order);
                    temp_tab_count = 0;
                    ++HOG_count;
                    temp_string = "";
                    i_storage.clear();
                }else if(temp_tab_count >= 3){ // first 3 tabs of OrthoFinder output do not contain prod_acc data but list the hog number, orthogroup, and parent clade
                    if(c == ' '){  // spaces are used after a comma. comma is used to delimit multiple genes from single sample in single HOG.
                        continue;
                    }
                    temp_string += c;
                }

            }
        }

        // Below is the diag print statements. 
        if(diag > 2){
            for(int i = 0; i < p; ++i){
                std::cout << "Feature table name" << '\t' << (*feature_tables_info)[i].name << '\t';
                std::cout << "Feature table CDS lines (no without protein)" << '\t' << (*feature_tables_info)[i].HOG.size() << '\n';
                if(diag > 2){
                    std::cout << "Gene_number (CDS & no without protein lines)" << '\t' << "HOG" << '\t' << "Product accession" << '\n';
                    for(int j = 0; j < (*feature_tables_info)[i].HOG.size(); ++j){
                        std::cout << j << '\t' << (*feature_tables_info)[i].HOG[j] << '\t' << (*feature_tables_info)[i].prod_acc[j] << '\t';
                        int match = (*feature_tables_info)[i].HOG[j];
                        if(match == -999){  // -999 was used to fill arrays/ vectors so we could recognize no data positions. 
                            std::cout << '\n';
                            continue;
                        }
                        for(int k = 0; k < p; ++k){
                            if(i == k){
                                continue;
                            }
                            for(int l = 0; l < (*feature_tables_info)[k].HOG.size(); ++l){
                                if((*feature_tables_info)[k].HOG[l] == match){
                                    std::cout << '\t' << (*feature_tables_info)[k].prod_acc[l] << '\t';
                                }
                            }
                            std::cout << '\n';
                        }

                    }
                    std::cout << '\n';
                }
            }
            if(diag > 2){
                std::cout << '\t' << '\t';
                for(int i = 0; i < p; ++i){  // p is the number of samples from the user provided input list
                    std::cout << (*feature_tables_info)[i].name << '\t' << "|" << '\t';
                }
                std::cout << '\n';

                for(int h = 0; h < HOG_count; ++h){
                    std::cout << "HOG" << '\t' << h << '\t';
                    for(int i = 0; i < p; ++i){ //p is total count of input files; i.e. how many feature tables
                        for(int j = 0; j < (*feature_tables_info)[i].HOG.size(); ++j){ 
                            if((*feature_tables_info)[i].HOG[j] == h){
                                std::cout << j << '\t' << (*feature_tables_info)[i].prod_acc[j] << '\t';
                            }
                        }
                        std::cout << "|" << '\t';
                    }
                    std::cout << '\n';
                }
            }

            if(diag > 2){

                std::vector<std::vector<std::string>> to_print_hog(p); // storage to print at end of this loop; p is number of user provided samples
                for(int i = 0; i < p; ++i){
                    to_print_hog[i].resize(p);
                    for(int j = 0; j < p; ++j){
                        to_print_hog[i][j] = "";
                    }
                }

                // Following calculates pairwise combinations. Could be moved to function & array "combo" passed to prevent repeat looping each time function called
                // n = 4; number of user provided samples. In this example the user provided 4. 
                // pairwise = (n(n-1))/2
                // pairwise = (4(3))/2 = 6 pairwise combinations in this example.
                // The combinations are:
                // 0    0-1
                // 1    0-2
                // 2    0-3
                // 3    1-2
                // 4    1-3
                // 5    2-3
                int samp_cnt = p;          // sample_count; number of user provided samples; better name
                int samp_trk = 1;          // sample tracker
                int loop_trk = 0;          // loop tracker
                int firt_com = 0;          // first combo. tracks "first" sample of combo
                int secd_com = 1;          // second combo. tracks "second" sample of combo
                int combo[pairwise][2];
                for(int i = 0; i < pairwise; ++i, ++loop_trk){
                    if(loop_trk == samp_cnt - samp_trk){
                        firt_com = samp_trk;
                        secd_com = firt_com + 1;
                        loop_trk = 0;
                        ++samp_trk;
                    }

                    combo[i][0] = firt_com;
                    combo[i][1] = secd_com;
                    ++secd_com;
                }

                for(int i = 0; i < pairwise; ++i){
                    std::cout << "COMBO" << '\t' << combo[i][0] << '\t' << combo[i][1] << '\n';
                }
                for(int i = 0; i < p; ++i){  // p is the number of samples from the user provided input list
                    std::cout << (*feature_tables_info)[i].name << '\t' << "|" << '\t';
                }
                std::cout << '\n';

                for(int b = 0; b < HOG_count; ++b){
                    if(diag > 0){
                        std::cout << "HOG" << '\t' << b << '\t';
                    }

                    for(int a = 0; a < pairwise; ++a){
                        for(int c = 0; c < 2; ++c){
                            for(int d = 0; d < hog_master[a][b][c].size(); ++d){
                                int e = hog_master[a][b][c][d];
                                if(e == -9999){
                                    continue;
                                }
                                bool already_added = 0; // make sure gene isn't already loaded to print
                                int home = combo[a][c]; // home will cycle through the index number of the user provided sample
                                std::string hog_string = std::to_string(e);
                                for(int i = 0; i < to_print_hog[home].size(); ++i){
                                    size_t found = to_print_hog[home][i].find(hog_string);
                                    if(found != std::string::npos){
                                        already_added = 1;
                                    }
                                }
                                if(already_added == 0){
                                    if(to_print_hog[home][0] == ""){
                                        to_print_hog[home][0] = hog_string + '\t' + (*feature_tables_info)[home].prod_acc[e] + '\t';
                                    }else{
                                        to_print_hog[home].push_back(hog_string + '\t' + (*feature_tables_info)[home].prod_acc[e] + '\t');
                                    }
                                }
                            }
                        }
                    }

                    if(diag > 0){
                        for(int i = 0; i < p; ++i){
                            for(int j = 0; j < to_print_hog[i].size(); ++j){
                                std::cout << to_print_hog[i][j];
                            }
                            std::cout << "|" << '\t';
                        }
                        std::cout << '\n';
                    }

                    for(int i = 0; i < p; ++i){ //reset contents of printing array
                        for(int j = 0; j < to_print_hog[i].size(); ++j){
                            to_print_hog[i][j] = "";
                        }
                    }
                }
            }
        }
        return hog_master;
    }
}

void assign_HOG(int temp_tab_count, int HOG_count, std::string temp_string, std::vector<strct_feature_tables_info> *feature_tables_info, std::vector<int> &i_storage, int benchmark, std::vector<std::string> name_storage, int ortho_finder_order){

    // Uncomment below statements for trouble shooting
    //std::cout << "TAB COUNT" << '\t' << temp_tab_count << '\n';
    //std::cout << "TEMP_STRING" << '\t' << temp_string << '\n';
    //std::cout << "FT FILE" << '\t' << (*feature_tables_info)[h].name << '\n';
    //std::cout << "FT SIZE ARRAY " << '\t' << (*feature_tables_info)[h].name << '\n';

    int h = temp_tab_count - 3; // -3 because the first 3 columns of OrthoFinder's output does not contain prod_acc data.
    if(ortho_finder_order == 0){
        for(int i = 0; i < (*feature_tables_info)[0].size_strct; ++i){
            std::string temp_ft_name = (*feature_tables_info)[i].name.substr(0, (*feature_tables_info)[i].name.find('_', 4)); // cut so only GCF_*. 
            std::string temp_N0_name = name_storage[h].substr(0, name_storage[h].find('_', 4));
            //std::cout << temp_ft_name << '\t' << temp_N0_name << '\n';
            if(temp_ft_name == temp_N0_name){
                h = i;
                goto stop_search;
            }
        }
        std::cout << "Error: user generated input file cannot find a match to GCF listed in N0.tsv. Let me know on Github so I can fix it or rerun with option ortho_finder_order -1" << '\n';
        std::exit(EXIT_FAILURE);
    }
    stop_search:;
    bool flag_i{0};             // Flag to warn that a prod_acc from the fasta file provided to OrthoFinder is not found in the feature table file.
    std::vector<int> i_stor;
    //need to change between .locus_tag and .prod_acc depending on input protein using one or the other
    if(benchmark == 1){
        for(int i = 0; i < (*feature_tables_info)[h].locus_tag.size(); ++i){
            if(temp_string == (*feature_tables_info)[h].locus_tag[i]){  // Match temp_string (prod_acc from OrthoFinder) to prod_acc from feature table
                flag_i = 1;
                (*feature_tables_info)[h].HOG[i] = HOG_count;  // Store the HOG into the feature table struct
                i_storage.push_back(i);  // Store the location / order from the feature table of the matching prod_acc
            }else{
                //std::cout << temp_string << '\t' << (*feature_tables_info)[h].locus_tag[i] << '\n';
            }
        }
    }else{ // benchmark == 0
        for(int i = 0; i < (*feature_tables_info)[h].prod_acc.size(); ++i){
            if(temp_string == (*feature_tables_info)[h].prod_acc[i]){  // Match temp_string (prod_acc from OrthoFinder) to prod_acc from feature table
                flag_i = 1;
                (*feature_tables_info)[h].HOG[i] = HOG_count;  // Store the HOG into the feature table struct
                i_storage.push_back(i);  // Store the location / order from the feature table of the matching prod_acc
            }else{
                //std::cout << temp_string << '\t' << (*feature_tables_info)[h].locus_tag[i] << '\n';
            }
        }
    }
    if(flag_i == 0){
        std::cout << "Warning: prod_acc " << temp_string << " from HOG file not found in " << (*feature_tables_info)[h].name << " feature table" << '\n';
    }
}  

void assign_HOG_match(int numb_samp, int temp_tab_count, int HOG_count, int pairwise, std::vector<int> i_storage, std::vector<int>*** hog_master, std::vector<strct_feature_tables_info> *feature_tables_info, std::vector<std::string> name_storage, int ortho_finder_order){

    // Following calculates pairwise combinations. Could be moved outside and passed in to prevent repeat looping each time function called
    // n = 4; number of user provided samples. In this example the user provided 4. 
    // pairwise = (n(n-1))/2
    // pairwise = (4(3))/2 = 6 pairwise combinations in this example.
    // The combinations are:
    // 0    0-1
    // 1    0-2
    // 2    0-3
    // 3    1-2
    // 4    1-3
    // 5    2-3
    int samp_cnt = numb_samp; // sample_count; number of user provided samples; better name
    int samp_trk = 1; // sample tracker
    int loop_trk = 0; // loop tracker
    int firt_com = 0; // first combo. tracks "first" sample of combo
    int secd_com = 1; // second combo. tracks "second" sample of combo
    int combo[pairwise][2];
    for(int i = 0; i < pairwise; ++i, ++loop_trk){
        if(loop_trk == samp_cnt - samp_trk){
            firt_com = samp_trk;
            secd_com = firt_com + 1;
            loop_trk = 0;
            ++samp_trk;
        }

        combo[i][0] = firt_com;
        combo[i][1] = secd_com;
        ++secd_com;
    }

    int home = temp_tab_count - 3;  // Which column from N0.tsv 
    if(ortho_finder_order == 0){
        for(int i = 0; i < (*feature_tables_info)[0].size_strct; ++i){
            std::string temp_ft_name = (*feature_tables_info)[i].name.substr(0, (*feature_tables_info)[i].name.find('_', 4)); // cut so only GCF_*. 
            std::string temp_N0_name = name_storage[home].substr(0, name_storage[home].find('_', 4));
            //std::cout << temp_ft_name << '\t' << temp_N0_name << '\n';
            if(temp_ft_name == temp_N0_name){
                home = i;
                goto stop_search;
            }
        }
        std::cout << "Error: user generated input file cannot find a match to GCF listed in N0.tsv. Let me know on Github so I can fix it or rerun with option ortho_finder_order -1" << '\n';
        std::exit(EXIT_FAILURE);
    }
    stop_search:;
    // Trouble shooting print statement
    /*
    for(int i = 0; i < pairwise; ++i, ++loop_trk){
        std::cout << combo[i][0] << '\t' << combo[i][1] << '\n';
    }
        
    if(hog_master[home][HOG_count][0][0] != -9999){
        std::cout << "Troubleshoot" << '\t' << home << '\t' << HOG_count << '\t' << hog_master[home][HOG_count][0][0] << '\n';
    }
    */

    // Loop assigns the order (location) of a prod_acc (gene) from the feature tables into hog_master
    // a is pairwise
    // HOG_count is the hog
    // b is tracking which "blast direction" or which sample of A vs. B to store in. 
    // b will allow to only compare sample A to B and keep a comparison of A vs. A from occuring
    // last dimension is if there mutiple matches on a single line of hog from input file

    if(i_storage.size() > 1){
        for(int d = 0; d < i_storage.size() - 1; ++d){
            for(int e = d + 1; e < i_storage.size(); ++e){
                if(i_storage[d] > i_storage[e]){  // changes order to be smallest to largest 
                    int temp = i_storage[e];
                    i_storage[e] = i_storage[d];
                    i_storage[d] = temp;
                }
            }
        }
    }
    
    for(int d : i_storage){ // store hog number
        for(int a = 0; a < pairwise; ++a){
            if(combo[a][0] == home){
                if(hog_master[a][HOG_count][0][0] == -9999){
                    hog_master[a][HOG_count][0][0] = d;
                }else{
                    hog_master[a][HOG_count][0].push_back(d);
                }
            }
            if(combo[a][1] == home){
                if(hog_master[a][HOG_count][1][0] == -9999){
                    hog_master[a][HOG_count][1][0] = d;
                }else{
                    hog_master[a][HOG_count][1].push_back(d);
                }
            }
        }
    }
}

std::vector<std::vector<std::string>> match_hog(std::vector<int>*** HOG_master, std::vector<strct_feature_tables_info> *feature_tables_info, int pairwise, int HOG_line, int window_size, double synteny_ratio, int diag, int ft_gene_together, int print_all_orthofinder, int run_all_orthofinder, int benchmark, int run_single_HOG, int prod_acc_flag){
    
    std::string out_file_name;
    out_file_name = "benchmark_submit_" + std::to_string(window_size) + '_' + std::to_string(synteny_ratio) + ".tsv";
    std::ofstream benchmark_file(out_file_name);
    std::vector<std::vector<std::string>> hog_match;
    hog_match.resize(HOG_line);
    for(int i = 0; i < HOG_line; ++i){
        hog_match[i].resize(1);
        for(int j = 0; j < 1; ++j){
            hog_match[i][j] = "";
        }
    }

    int samp_cnt = (*feature_tables_info)[0].size_strct;  // sample count
    int samp_trk = 1; // sample tracker
    int loop_trk = 0; // loop tracker
    int firt_com = 0; // first combo. tracks "first" sample of combo
    int secd_com = 1; // second combo. tracks "second" sample of combo
    int combo[pairwise][2];
    for(int i = 0; i < pairwise; ++i, ++loop_trk){
        if(loop_trk == samp_cnt - samp_trk){
            firt_com = samp_trk;
            secd_com = firt_com + 1;
            loop_trk = 0;
            ++samp_trk;
        }

        combo[i][0] = firt_com;
        combo[i][1] = secd_com;
        ++secd_com;
    }

    /*
    std::cout << "COMBO" << '\n';
    for(int i = 0; i < pairwise; ++i){
        std::cout << combo[i][0] << '\t' << combo[i][1] << '\n';
    }
    */


    // I and J were swapped to fix a bug (needed to print results from pairwise match together)
    for(int j = 0; j < HOG_line; ++j){
        if(run_single_HOG != -1){
            if(run_single_HOG > HOG_line || run_single_HOG < 0){
                std::cout << "User requested to run a single HOG outside the number of HOGS in OrthoFinder file!" << '\n';
                std::exit(EXIT_FAILURE);
            }
            if(j != run_single_HOG){
                continue;
            }
        }
        int sum = 0;
        if(HOG_master[0][j][0][0] != -9999){
            sum += HOG_master[0][j][0].size();
        }
        if(HOG_master[0][j][1][0] != -9999){
            sum += HOG_master[0][j][1].size();
        }
        for(int i = 0; i < samp_cnt - 2; ++i){
            if(HOG_master[combo[i][1]][j][1][0] != -9999){
                sum += HOG_master[combo[i][1]][j][1].size();
            }
           
        }
        if(diag > 3){
            std::cout << "SUM" << '\t' << j << '\t' << sum << '\n';
        }

        std::vector<std::vector<std::string>> locus_tag(500, std::vector<std::string>(2)); //10000
        std::vector<std::vector<int>> locus_ft(500, std::vector<int>(2));  // which feature table (index) does the locus tag belong to
        std::vector<std::vector<std::string>> prod_acc(500, std::vector<std::string>(2));

        bool run_alg_flag{0};  // If any genome member of the HOG has more than 2 genes, then run the alg on all members of the HOG for print statements

        for(int i = 0; i < pairwise; ++i, ++loop_trk){
            if(run_all_orthofinder == 0){
                if(HOG_master[i][j][0].size() > 1 || HOG_master[i][j][1].size() > 1){
                    run_alg_flag = 1;
                }
            }else{
                run_alg_flag = 1;
                break;
            }
        }
        if(run_alg_flag == 1){  
            int vec_loop_cnt = 0; // vector loop counter to know what index to store the vector, temp_prod_vector, into prod_acc
            for(int i = 0; i < pairwise; ++i){
            
                if(diag > 0){
                    std::cout << "STARTING HOG" << '\t' << j << std::endl;
                }

                int bigger{-9999}; // stores number of entries of the larger 2 from hog_master[i][j][0] vs hog_master[i][j][1]
                int countdown{-9999};  // max number of genes contributed from one pairwise sample at a particular hog
                int hmk{-9999};  // hog_master k
                int hmnk{-9999}; // hog_master not k
                int ftk{-9999};  // feature table k
                int ftnk{-9999}; // feature table not k
                if(HOG_master[i][j][0].size() >= HOG_master[i][j][1].size()){
                    bigger = HOG_master[i][j][0].size();  // 0
                    countdown = HOG_master[i][j][1].size();  // 1
                    ftk = combo[i][1]; // 1
                    ftnk = combo[i][0]; // 0
                    hmk = 1;  // 1
                    hmnk = 0;  // 0 
                }else{
                    bigger = HOG_master[i][j][1].size();  // 1
                    countdown = HOG_master[i][j][0].size();  // 0             
                    ftk = combo[i][0]; // 0
                    ftnk = combo[i][1]; // 1
                    hmk = 0;  // 0
                    hmnk = 1;  // 1
                }

                // The together vectors are used when the -together option is set to 1 (the default). See --help on the -together option
                std::vector<std::vector<int>> together1(countdown, std::vector<int> (1, -9999));
                std::vector<std::vector<int>> together2(HOG_master[i][j][hmnk].size(), std::vector<int> (1, -9999));
                if(ft_gene_together == 1){  // -togther is set to 1, the default
                    int tg_count{1}; // together counter
                    for(int l = 0; l < countdown; ++l){
                        int m = HOG_master[i][j][hmk][l];
                        if(m == -9999){
                            continue;
                        } 
                        for(int q = 0; q < countdown; ++q){
                            if(q == l){
                                continue;
                            }
                            int m2 = HOG_master[i][j][hmk][q];
                            if(m2 == -9999){
                                continue;
                            }

                            if(abs(m - m2) <= tg_count){ // checks if 2 genes in the same HOG are sequential in first genome

                                ++tg_count;
                                if(diag > 2){
                                    std::cout << "STORING TOGETHER 1" << '\t' << i << '\t' << j << '\t' << hmk << '\t' << l << '\t' << q << '\n';
                                    std::cout << "FT1 INFO" << '\t' << (*feature_tables_info)[ftk].locus_tag[m] << std::endl;
                                }
                                if(together1[l][0] == -9999){
                                    together1[l][0] = l;
                                }
                                together1[l].push_back(q);

                            }
                        }
                    }
                    for(int l = 0; l < countdown; ++l){
                        for(int q = 0; q < together1[l].size() - 1; ++q){
                            if(together1[l][q] > together1[l][q + 1]){ // checks if genes were stored in order
                                if(diag > 2){
                                    std::cout << "REORDER TOGETHER1" << '\n';
                                }
                                int temp = together1[l][q + 1];
                                together1[l][q + 1] = together1[l][q];
                                together1[l][q] = temp;
                            }else if (diag > 2){
                                std::cout << "NO REORDER TOGETHER" << '\t' << together1[l][q] << '\t' << together2[l][q + 1] << '\n';
                            }
                        }
                    }
                    tg_count = 1;
                    for(int p = 0; p < HOG_master[i][j][hmnk].size(); ++p){
                        int t1 = HOG_master[i][j][hmnk][p];  // test1;
                        if(t1 == -9999){
                            continue;
                        }
                        for(int q = 0; q < HOG_master[i][j][hmnk].size(); ++q){
                            int t2 = HOG_master[i][j][hmnk][q];  // test2
                            if(t2 == -9999 || q == p){
                                continue;
                            }
                            if(diag > 2){
                                std::cout << "Chceking if together 2" << '\t' << p << '\t' << t1 << '\t' << (*feature_tables_info)[ftnk].locus_tag[t1] << '\t' << q << '\t' << t2 << '\t' << (*feature_tables_info)[ftnk].locus_tag[t2] << '\t' << together2[p][0] << '\n';
                            }
                            if(abs(t1 - t2) <= tg_count){ // checks if 2 genes in the same HOG are sequential in other genome (second)
                                ++tg_count;
                                if(diag > 2){
                                    std::cout << "STORING TOGETHER 2" << '\t' << "temp" << '\t' << i << '\t' << j << '\t' << hmnk << '\t' <<p << '\t' << q << '\n';
                                    std::cout << "FT2 INFO" << '\t' << (*feature_tables_info)[ftnk].locus_tag[t1] << std::endl;
                                }
                                if(together2[p][0] == -9999){
                                    together2[p][0] = p;
                                }
                                together2[p].push_back(q);
                            }
                            if(diag > 2){
                                std::cout << "Together 2" << '\t' << p << '\t' << q << '\t' << together2[p][0] << '\n';
                            }
                        }
                    }
                    for(int l = 0; l < HOG_master[i][j][hmnk].size(); ++l){
                        for(int q = 0; q < together2[l].size() - 1; ++q){
                            if(together2[l][q] > together2[l][q + 1]){ // checks if genes were stored in order
                                if(diag > 2){
                                    std::cout << "SWITCH TOGETHER2" << '\n';
                                }
                                int temp = together2[l][q + 1];
                                together2[l][q + 1] = together2[l][q];
                                together2[l][q] = temp;
                            }else if (diag > 2){
                                std::cout << "NO SWITCH TOGETHER 2" << '\t' << together2[l][q] << '\t' << together2[l][q+1] << '\n';
                            }
                        }
                    }
                }
                if(diag > 2){
                    std::cout << "CHECKING TOGETH 1 SORT" << '\n';
                    for(int l = 0 ; l < countdown; ++l){
                        for(int q = 0 ; q < together1[l].size(); ++q){
                            std::cout << together1[l][q] << '\t';
                        }
                        std::cout << '\n';
                    }
                    std::cout << "CHECKING TOGETHER 2 SORT" << '\n';
                    for(int l = 0; l < HOG_master[i][j][hmnk].size(); ++l){
                        for(int q = 0; q < together2[l].size(); ++q){
                            std::cout << together2[l][q] << '\t';
                        }
                        std::cout << '\n';
                    }
                    std::cout << "FINISH CHECK TOGETHER SORT" << '\n';
                }

                for(int l = 0; l < countdown; ++l){
                    if(diag > 1){
                        std::cout << "Starting countdown" << '\t' << "i" << '\t' << i << '\t' << "j" << '\t' << j << '\t' << "L" << '\t' << l << std::endl;
                    }
                    int m = HOG_master[i][j][hmk][l];
                    if(m == -9999){
                        continue;
                    }

                    int **zstore = new int*[window_size + 1];
                    int **ustore = new int*[window_size + 1]; // {-999}
                    int zcount[HOG_master[i][j][hmnk].size()]{0};
                    int z_neg_count[HOG_master[i][j][hmnk].size()]{0};  // Tracks how many -999 
                    int u_neg_count[HOG_master[i][j][hmnk].size()]{0};
                    int z_scaff_count[HOG_master[i][j][hmnk].size()]{0}; // Tracks scaffold count
                    int u_scaff_count[HOG_master[i][j][hmnk].size()]{0};
                    for(int p = 0; p < HOG_master[i][j][hmnk].size(); ++p){

                        int s = HOG_master[i][j][hmnk][p];
                        if(s == -9999){
                            continue;
                        }

                        int z_max_operon_count{1};
                        int u_max_operon_count{1};
                        while_loop(HOG_master, j, p, zstore, z_neg_count, z_scaff_count, feature_tables_info, ftk, m, window_size, diag, ftnk, pairwise, combo, z_max_operon_count); // Loops through each window of the 2 genomes to identify what genes to compare
                        while_loop(HOG_master, j, p, ustore, u_neg_count, u_scaff_count, feature_tables_info, ftnk, s, window_size, diag, ftk, pairwise, combo, u_max_operon_count);

                        if(diag > 1){
                            std::cout << '\n';
                        }
                        bool ucheck[window_size + 1]{0};  // prevents the same u from matching mutiple z
                        bool zcheck[window_size + 1]{0};  // prevents the same u from matching mutiple z
                        for(int z = 0; z < window_size + 1; ++z){ // counts how many matches were in a window
                            for(int u = 0; u < window_size + 1; ++u){
                                for(int z2D = 0; z2D < z_max_operon_count; ++z2D){
                                    if(zstore[z][z2D] == -999){
                                        continue;
                                    }
                                    for(int u2D = 0; u2D < u_max_operon_count; ++u2D){
                                        if(ustore[u][u2D] == -999){
                                            continue;
                                        }
                                        if(zstore[z][z2D] == ustore[u][u2D] && ucheck[u] != 1 && zcheck[z] != 1 && ustore[u][u2D] != -999 && zstore[z][z2D] != -999){ 
                                            if(diag > 0){
                                                std::cout << "Match" << '\t' << zstore[z][z2D] << '\t' << ustore[u][u2D] << '\t' << zcheck[z] << '\t' << ucheck[u] << '\t' << z << '\t' << u << std::endl;
                                            }
                                            zcount[p] += 1;
                                            ucheck[u] = 1;
                                            zcheck[z] = 1;
                                            break;
                                        }else{
                                            //std::cout << zstore[z][z2D] << '\t' << ustore[u][u2D] << '\t' << "ucheck" << '\t' << ucheck[u] << '\n';
                                        }
                                    }
                                }
                            }
                        }

                        if(diag > 0){
                            std::cout << "Final matching count:" << '\t' << zcount[p] << std::endl;
                        }
                    }

                    double z_max = synteny_ratio;
                    std::vector<int> z_ind(1, -9999);  // z index
                    std::vector<int> p_ortho(1, -9999); // print orthofinder, for benchmarking option
                    for(int p = 0; p < HOG_master[i][j][hmnk].size(); ++p){
                        int reduce_ratio{0}; // If zstore or ustore contains -999 (from a scaffold edge being in the window), reduce the denominator of the required ratio below
                        if(z_scaff_count[p] > u_scaff_count[p]){
                            reduce_ratio = z_scaff_count[p];
                        }else{
                            reduce_ratio = u_scaff_count[p];
                        }
                        int adjst_window = window_size - reduce_ratio;
                        if(adjst_window == 0){  // If there was a scaffold immediately on both sides of a gene allowing for no comparison
                            continue;
                        }
                        double z_ratio = double(zcount[p]) / adjst_window;
                        if(diag > 2){
                            std::cout << "Z count" << '\t' << zcount[p] << '\t' << "adjust window" << '\t' << adjst_window << '\n';
                        }

                        if(z_ratio > z_max && zcount[p] > 0){ // if the z_ratio is greater than a previous max
                            bool erased{0};
                            if(z_ind.size() > 1 && print_all_orthofinder == 0){  // remove extra elements that matched previous smaller size if we are not printing all orthofinder
                                z_ind.erase(z_ind.begin() + 1, z_ind.end());
                                erased = 1;
                            }
                            if(erased == 0 && z_ind.size() > 1){
                                z_ind.push_back(p);
                            }else{
                                z_ind[0] = p;
                            }

                            if(diag > 0){
                                int temp = HOG_master[i][j][hmnk][z_ind[z_ind.size() - 1]];
                                std::cout << "SOG BETTER" << '\t' << "HOG" << '\t' << j << '\t' << (*feature_tables_info)[ftk].locus_tag[m] << '\t' << (*feature_tables_info)[ftnk].locus_tag[temp]<<  '\t' << "ratio" << '\t' << z_ratio << '\t' << "max" << '\t' << z_max << '\t' << "adjust window" << '\t' <<   adjst_window <<  '\t' << "reduce ratio" << '\t'<<    reduce_ratio << '\n';
                            }

                            z_max = z_ratio; 
                        }else if((abs(z_ratio - z_max)) < 0.0001){ // if z ratio matches previous z max
                            
                            if(z_ind[0] == -9999){
                                z_ind[0] = p;
                            }else{
                                z_ind.push_back(p);
                            }
                            if(diag > 0){
                                int temp = HOG_master[i][j][hmnk][z_ind[z_ind.size() - 1]];
                                std::cout << "Z SAME, PUSH BACK" << '\t' << "p" << '\t' << p << '\t' << "zcount[p]" << '\t' << zcount[p] << '\t' << "zmax" << '\t' << z_max<< '\t' << "z_ratio" << '\t' << z_ratio << '\t' << (*feature_tables_info)[ftk].locus_tag[m] << '\t' << (*feature_tables_info)[ftnk].locus_tag[temp] <<  '\n';
                            }
                        }else if(print_all_orthofinder == 1 && HOG_master[i][j][hmnk][p] != -9999){ // if we are printing all orthofinder (printing even if no improvement could be made to HOG)
                            if(benchmark == 0){ // prevents z_ind being inserted when benchmarking (benchmark == 1) as this makes the file have incorrect orthologs. cannot run outfile and benchmark at same time
                                if(z_ind[0] == -9999){
                                    z_ind[0] = p;
                                }else{
                                    z_ind.push_back(p);
                                }
                            }
                            if(p_ortho[0] == -9999){
                                p_ortho[0] = p;
                            }else{
                                p_ortho.push_back(p);
                            }
                        }
                        if(z_ratio < z_max && diag > 0 && HOG_master[i][j][hmnk][p] != -9999){
                            int temp = HOG_master[i][j][hmnk][p];
                            std::cout << "SOG NOT BETTER" << '\t' << "HOG" << '\t' << j << '\t' << (*feature_tables_info)[ftk].locus_tag[m] << '\t' << (*feature_tables_info)[ftnk].locus_tag[temp]<<  '\t' << "ratio" << '\t' << z_ratio << '\t' << "max" << '\t' << z_max << '\t' << "adjust window" << '\t' <<   adjst_window <<  '\t' << "reduce ratio" << '\t'<<    reduce_ratio << '\t' << "FTK" << '\t' << ftk << '\t' << "FTNK" << '\t' << ftnk << '\n';
                        }
                    }

                    if(p_ortho[0] != -9999 && z_ind[0] == -9999){ // If we should print benchmark but there was no synteny support and we want to print all orthofinder
                        std::string ft1_name_locus_tag = (*feature_tables_info)[ftk].locus_tag[m];
                        if(benchmark == 1){
                            for(int q = 0; q < p_ortho.size(); ++q){
                                int temp = HOG_master[i][j][hmnk][p_ortho[q]];
                                std::string ft2_name_locus_tag = (*feature_tables_info)[ftnk].locus_tag[temp];
                                if(diag > 2)std::cout << "Storing into benchmark file" << '\t' << ft1_name_locus_tag << '\t' << ft2_name_locus_tag << '\n';
                                benchmark_file << ft1_name_locus_tag << '\t' << ft2_name_locus_tag << '\n';
                            }
                        }
                        continue;
                    }else if(z_ind[0] == -9999){
                        continue;
                    }
                    std::string ft1_name_locus_tag = (*feature_tables_info)[ftk].locus_tag[m];
                    std::string ft1_name_prod_acc = (*feature_tables_info)[ftk].prod_acc[m];
                    if(together1[l][0] != -9999){ // Stores two sequential genes from same HOG together
                        ft1_name_locus_tag = "";
                        ft1_name_prod_acc = "";
                        for(int q = 0; q < together1[l].size(); ++q){
                            int qtemp = together1[l][q];
                            int mtemp = HOG_master[i][j][hmk][qtemp];
                            if(ft1_name_locus_tag != ""){
                                ft1_name_locus_tag += '\t';
                                ft1_name_prod_acc += '\t';
                            }
                            ft1_name_locus_tag += (*feature_tables_info)[ftk].locus_tag[mtemp];
                            ft1_name_prod_acc += (*feature_tables_info)[ftk].prod_acc[mtemp];
                        }
                    }
                    if(diag > 1){
                        std::cout << "ft1 name locus tag" << '\t' << ft1_name_locus_tag << '\n';
                    }
                    int larger_zind_or_together2{0};
                    if(z_ind.size() > together2.size()){
                        larger_zind_or_together2 = z_ind.size();
                    }else{
                        larger_zind_or_together2 = together2.size();
                    }
                    std::string verify_below[larger_zind_or_together2]; // verify the locus tag wasn't already stored to print in loop above from being sequential in the genome. Will happen if 2 genes have same best score and are sequential in order
                    for(int q = 0; q < larger_zind_or_together2; ++q){
                        verify_below[q] = "";
                    }
                    if(diag > 1){
                        std::cout << "value stored at z_ind[0]" << '\t' << z_ind[0] << '\t' << "size of z_ind" << '\t' << z_ind.size() << '\n';
                        std::cout << "values stored in z_ind are" << '\n';
                        for(int hgh = 0; hgh < z_ind.size(); ++hgh){
                            std::cout << z_ind[hgh] << '\n';
                        }
                        std::cout << "end values stored in z_ind" << '\n';
                    }

                    int temp = HOG_master[i][j][hmnk][z_ind[0]];
                    std::string ft2_name_locus_tag = (*feature_tables_info)[ftnk].locus_tag[temp];
                    std::string ft2_name_prod_acc = (*feature_tables_info)[ftnk].prod_acc[temp];
                    verify_below[0] = (*feature_tables_info)[ftnk].locus_tag[temp];
                    if(diag > 1){
                        std::cout << "ft2 name locus tag" << '\t' << ft2_name_locus_tag << '\n'; // ft_name_locus_tag
                    }
                    int qtemp_store{-1}; 
                    for(int q = 0; q < together2[z_ind[0]].size(); ++q){ // if a "first" gene was stored in together 2, we need to know where (index) to store it below 
                        if(HOG_master[i][j][hmnk][together2[z_ind[0]][q]] == temp && together2[z_ind[0]][q] != -9999){ // loop through together 2
                            qtemp_store = together2[z_ind[0]][q]; 
                            if(diag > 2)std::cout << "breaking on qtemp match" << '\t' << HOG_master[i][j][hmnk][together2[z_ind[0]][q]] << '\t' << temp << '\n';
                            break;
                        }else{
                            //std::cout << HOG_master[i][j][hmnk][q] << '\t' << temp << '\t' << together2[z_ind[0]][q] << '\t' << HOG_master[i][j][hmnk][together2[z_ind[0]][q]] << '\n';
                        }
                    }
                    if(diag > 2)std::cout << "qtemp_store" << '\t' << qtemp_store << '\n';
                    std::vector<std::string> ft_2_temp_locus(HOG_master[i][j][hmnk].size());
                    std::vector<std::string> ft_2_temp_prod(HOG_master[i][j][hmnk].size());
                    if(qtemp_store == -1){ // no gene was stored in together 2
                        ft_2_temp_locus[0] = ft2_name_locus_tag; // store "first" gene at index 0
                        ft_2_temp_prod[0] = ft2_name_prod_acc;
                    }else if(qtemp_store > HOG_master[i][j][hmnk].size()){ //safety check
                        std::cout << "qtempstore store greater than HOG_master[i][j][hmnk].size()" << '\n';
                        std::exit(EXIT_FAILURE);
                    }else{ // "first" gene was stored in together 2, store it at index qtemp_store
                        ft_2_temp_locus[qtemp_store] = ft2_name_locus_tag;
                        ft_2_temp_prod[qtemp_store] = ft2_name_prod_acc;
                    }
                    
                    for(int q = 0; q < z_ind.size(); ++q){  
                        if(together2[z_ind[q]][0] != -9999){  // check if another locus tag ("second" gene) was sequential to best match in feature_table; if so - also print it
                            if(diag > 1){
                                std::cout << "printing values in together2" << '\n';
                                for(int hgh = 0; hgh < together2[z_ind[q]].size(); ++hgh){
                                    std::cout << together2[z_ind[q]][hgh] << '\n';
                                }
                                std::cout << "finish printing values in together 2" << '\n';
                            }
                            ft2_name_locus_tag = "";
                            ft2_name_prod_acc = "";

                            for(int r = 0; r < together2[z_ind[q]].size(); ++r){
                                int qtemp = together2[z_ind[q]][r];                                
                                int mtemp = HOG_master[i][j][hmnk][qtemp];
                                std::string temp_temp = (*feature_tables_info)[ftnk].locus_tag[mtemp];
                                ft_2_temp_locus[qtemp] = temp_temp;
                                temp_temp = (*feature_tables_info)[ftnk].prod_acc[mtemp];
                                ft_2_temp_prod[qtemp] = temp_temp;
                                verify_below[qtemp] = (*feature_tables_info)[ftnk].locus_tag[mtemp];
                            }
                        }
                    }
                    ft2_name_locus_tag = "";
                    ft2_name_prod_acc = "";
                    bool first_entry{0};
                    for(int q = 0; q < ft_2_temp_locus.size(); ++q){
                        if(diag > 2){
                            std::cout << ft_2_temp_locus[q] << '\n';
                        }
                        if(ft_2_temp_locus[q] != ""){
                            if(first_entry == 0){ // only inserts a tab after the first entry
                                first_entry = 1;
                            }else{
                                ft2_name_locus_tag += '\t';
                                ft2_name_prod_acc += '\t';
                            }
                            ft2_name_locus_tag += ft_2_temp_locus[q];
                            ft2_name_prod_acc += ft_2_temp_prod[q];
                        }
                    }
                    if(diag > 1){
                        std::cout << "FT2 name locus tag post together" << '\t' << ft2_name_locus_tag << '\n';
                    }
                    if(z_ind.size() > 1){  // if best match in hog had same synteny_ratio as other possible matches within HOG, write all best possible
                        for(int q = 1; q < z_ind.size(); ++q){
                            bool check{0};
                            // verify the locus tag wasn't already stored to print in loop above from being sequential in the genome. Will happen if 2 genes have same best score and are sequential in order
                            temp = HOG_master[i][j][hmnk][z_ind[q]];
                            for(int r = 0; r < larger_zind_or_together2; ++r){
                                if(verify_below[r] == (*feature_tables_info)[ftnk].locus_tag[temp]){
                                    check = 1;
                                    if(diag > 1){
                                        std::cout << "breaking on requiring verify duplication" << '\n';
                                    }
                                    break;
                                }
                            }
                            if(check == 0 && diag > 1){
                                std::cout << "check is 0" << '\n';
                                for(int r = 0; r < larger_zind_or_together2; ++r){
                                    std::cout << "verify below" << '\t' << r << '\t' << verify_below[r] << '\t' << "ft_locus" << '\t' << (*feature_tables_info)[ftnk].locus_tag[temp] << '\n';
                                }
                            }else if(check == 1 && diag > 1){
                                std::cout << "check is 1" << '\n';
                                for(int r = 0; r < larger_zind_or_together2; ++r){
                                    std::cout << "verify below" << '\t' << r << '\t' << verify_below[r] << '\t' << "ft_locus" << '\t' << (*feature_tables_info)[ftnk].locus_tag[temp] << '\n';
                                }
                            }
                            if(check == 0){
                                ft2_name_locus_tag += '\t';
                                ft2_name_locus_tag += (*feature_tables_info)[ftnk].locus_tag[temp];
                                ft2_name_prod_acc += '\t';
                                ft2_name_prod_acc += (*feature_tables_info)[ftnk].prod_acc[temp];
                            }
                        }
                    }

                    std::vector<int> temp_locus_ft(2); 
                    std::vector<std::string> temp_locus_vector(2);
                    std::vector<std::string> temp_prod_vector(2);
                    if(hmk == 0){
                        temp_locus_vector[0] = ft1_name_locus_tag;
                        temp_locus_vector[1] = ft2_name_locus_tag;
                        temp_prod_vector[0] = ft1_name_prod_acc;
                        temp_prod_vector[1] = ft2_name_prod_acc;
                        temp_locus_ft[0] = ftk;
                        temp_locus_ft[1] = ftnk;
                    }else if(hmk == 1){
                        temp_locus_vector[0] = ft2_name_locus_tag;
                        temp_locus_vector[1] = ft1_name_locus_tag;
                        temp_prod_vector[0] = ft2_name_prod_acc;
                        temp_prod_vector[1] = ft1_name_prod_acc;
                        temp_locus_ft[0] = ftnk; 
                        temp_locus_ft[1] = ftk; 
                    }
                    if(diag > 3){
                        std::cout << "Printing temp locus tag and prod acc vectors" << '\t' << j << '\t' << i << '\t' << l << '\t' << temp_locus_vector[0] << '\t' << temp_locus_vector[1] << '\t' << temp_prod_vector[0] << '\t' << temp_prod_vector[1] << '\t' << '\n';
                        if(hmk == 0){
                            std::cout << "ftk" << '\t' << ftk << '\t' << "ftnk" << '\t' << ftnk << '\n';
                        }else{
                            std::cout << "ftnk" << '\t' << ftnk << '\t' << "ftk" << '\t' << ftk << '\n';
                        }
                    }

                    if(benchmark == 1){ // printing to benchmark file
                        std::vector<std::string> bench1;
                        std::vector<std::string> bench2;
                        std::string bench = "";
                        for(int xx = 0; xx < ft1_name_locus_tag.size(); ++xx){
                            if(ft1_name_locus_tag[xx] != '\t'){
                                bench += ft1_name_locus_tag[xx];
                            }else if(ft1_name_locus_tag[xx] == '\t'){
                                bench1.push_back(bench);
                                bench = "";
                            }
                            if(xx == ft1_name_locus_tag.size() - 1){
                                bench1.push_back(bench);
                                bench = "";
                            }
                        }
                        bench = "";
                        for(int xx = 0; xx < ft2_name_locus_tag.size(); ++xx){
                            if(ft2_name_locus_tag[xx] != '\t'){
                                bench += ft2_name_locus_tag[xx];
                            }else if(ft2_name_locus_tag[xx] == '\t'){
                                bench2.push_back(bench);
                                bench = "";
                            }
                            if(xx == ft2_name_locus_tag.size() - 1){
                                bench2.push_back(bench);
                                bench = "";
                            }
                        }
                        for(int xx = 0; xx < bench1.size(); ++xx){
                            for(int yy = 0 ; yy < bench2.size(); ++yy){
                                benchmark_file << bench1[xx] << '\t' << bench2[yy] << '\n';
                            }
                        }
                    }

                    if(vec_loop_cnt >= 500){ // safety check
                        std::cout << "vec loop count over 500" << '\n';
                        std::exit(EXIT_FAILURE);
                    }

                    for(int q = 0; q < vec_loop_cnt; ++q){ // verify pairwise orthologs were not already written to be printed
                        if(locus_tag[q][0] == temp_locus_vector[0]){
                            if(locus_tag[q][1] == temp_locus_vector[1]){
                                goto skip_vec_loop_count;
                            }
                        }else if(locus_tag[q][0] == temp_locus_vector[1]){
                            if(locus_tag[q][1] == temp_locus_vector [0]){
                                goto skip_vec_loop_count;
                            }
                        }else{

                        }
                    }

                    locus_tag[vec_loop_cnt] = temp_locus_vector;
                    locus_ft[vec_loop_cnt] = temp_locus_ft;
                    prod_acc[vec_loop_cnt]= temp_prod_vector;
                    ++vec_loop_cnt;

                    skip_vec_loop_count:;

                    if(diag > 0){
                        std::cout << "HOG" << '\t' << j << '\t' << "z_ind" << '\t' << z_ind[0] << '\t' << together2[z_ind[0]][0] << '\n';
                    }

                    for(int del = 0; del < window_size + 1; ++del){
                            delete[] zstore[del];
                            delete[] ustore[del];
                    }
                    delete[] zstore;
                    delete[] ustore;

                }

                if(diag > 0){
                    std::cout << "ENDING HOG" << '\t' << j << std::endl;
                }

            }

            if(diag > 1){
                std::cout << "locus_tag.size() pre lt_size count and adjustment" << '\t' << locus_tag.size() << '\n';
            }
            int lt_size{0};
            int pa_size{0};
            for(int a = 0; a < locus_tag.size(); ++a){
                if(locus_tag[a][0] != ""){
                    ++lt_size;
                }else{
                    break;
                }
            }
            for(int a = 0; a < prod_acc.size(); ++a){
                if(prod_acc[a][0] != ""){
                    ++pa_size;
                }else{
                    break;
                }
            }
            locus_tag.resize(lt_size);
            prod_acc.resize(pa_size);

            //DO NOT DELETE CODE COMMENTED OUT BELOW
            if(diag > 0){
                std::cout << "j" << '\t' << j << '\t' << "prod_acc_size" << '\t' << prod_acc.size() << '\t' << prod_acc[0].size() << '\t' << "lt_size" << '\t' <<  lt_size << '\n'; 
                std::cout << "PRINT OLD" << '\n';
                for(int xx = 0; xx < locus_tag.size(); ++xx){
                    for(int yy = 0; yy < locus_tag[xx].size(); ++yy){
                        std::cout << locus_tag[xx][yy] << '\t' << "$$$" << '\t' << locus_tag.size() << '\t' << locus_tag[xx].size() << '\t';
                    }
                std::cout << '\n' << '\n';
                }
            
                //for(int xx = 0; xx < prod_acc.size(); ++xx){
                //    for(int yy = 0; yy < prod_acc[xx].size(); ++yy){
                //        std::cout << prod_acc[xx][yy] << '\t' << "$$$" << '\t' << prod_acc.size() << '\t' << prod_acc[xx].size() << '\t';
                //    }
                //    std::cout << '\n';
                //}
            }

            if(diag > 1){
                std::cout << "starting match" << std::endl;
                for(int xx = 0; xx < locus_tag.size(); ++xx){
                    for(int yy = 0; yy < locus_tag[xx].size(); ++yy){
                        std::cout << locus_tag[xx][yy] << '\t' << "^^^" << '\t' << locus_tag.size() << '\t' << locus_tag[xx].size() << '\t';
                    }
                    std::cout << '\n' << '\n';
                }
            }

            auto start = std::chrono::high_resolution_clock::now();
            std::vector<std::vector<std::vector<std::string>>> final_answer(locus_tag.size(), std::vector<std::vector<std::string>> ((*feature_tables_info)[0].size_strct));            
            std::vector<std::vector<std::vector<std::string>>> locus_tag_apart(locus_tag.size(), std::vector<std::vector<std::string>> (2));
            if(prod_acc_flag == 1){
                for(int z = 0 ; z < locus_tag.size(); ++z){
                    for(int zz = 0; zz < 2; ++zz){
                        bool found_t_z{0};
                        size_t pos_z{0};
                        std::string search_z = "";
                        while(pos_z != std::string::npos){
                            search_z = "";
                            find_substr(found_t_z, pos_z, locus_tag, z, zz, search_z); // split locus tag from one element of many tabs to many elements with no tabs
                            locus_tag_apart[z][zz].push_back(search_z);
                        }
                    }
                }
            }else{
                for(int z = 0 ; z < prod_acc.size(); ++z){
                    for(int zz = 0; zz < 2; ++zz){
                        bool found_t_z{0};
                        size_t pos_z{0};
                        std::string search_z = "";
                        while(pos_z != std::string::npos){
                            search_z = "";
                            find_substr(found_t_z, pos_z, prod_acc, z, zz, search_z); // split locus tag from one element of many tabs to many elements with no tabs
                            locus_tag_apart[z][zz].push_back(search_z);
                        }
                    }
                }
            }

            if(diag > 1){
                std::cout << "LOCUS TAG APART" << '\n';
                for(int z = 0; z < locus_tag.size(); ++z){
                    for(int zz = 0; zz < 2; ++zz){
                        for(int zzz = 0; zzz < locus_tag_apart[z][zz].size(); ++zzz){
                            std::cout << locus_tag_apart[z][zz][zzz] << '\t';
                        }

                    }
                    std::cout << '\n';
                }
            }

            // this nested for loop combines ortholog pairs. If A matches B and A matches C, create single entry of A-B-C instead of 3 entries of A-B, B-C, A-C
            for(int a = 0; a < locus_tag_apart.size(); ++a){
                if(diag > 1)std::cout << "starting a" << '\t' << a << '\n';
                bool do_assign{1};
                for(int aa = 0; aa < locus_tag_apart[a].size(); ++aa){ // If A-B is already in final answer and A-C exist, insert C onto A-B
                    if(diag > 1)std::cout << "starting aa" << '\t' << aa << '\t' << locus_tag_apart[a][aa][0] << '\n';
                    for(int aaa = 0; aaa < locus_tag_apart[a][aa].size(); ++aaa){
                        if(locus_tag_apart[a][aa][aaa] == ""){
                            std::cout << "BLANK A AA AAA" << '\n'; // safety check, should never trigger
                        }
                        if(diag > 1)std::cout << "starting aaa" << '\t' << aaa << '\t' << locus_tag_apart[a][aa][aaa] <<'\n';
                        for(int x = 0; x < final_answer.size(); ++x){
                            for(int xx = 0; xx < final_answer[x].size(); ++xx){
                                for(int xxx = 0; xxx < final_answer[x][xx].size(); ++xxx){
                                    if(final_answer[x][xx][xxx] == ""){
                                        std::cout << "BLANK X XX XXX" << '\n';
                                    }
                                    if(locus_tag_apart[a][aa][aaa] == final_answer[x][xx][xxx]){
                                        if(diag > 1)std::cout << "FIRST MATCH" << '\t' << locus_tag_apart[a][aa][aaa] << '\n';
                                        for(int bb = 0; bb < locus_tag_apart[a].size(); ++bb){
                                            int pos = locus_ft[a][bb];
                                            for(int bbb = 0; bbb < locus_tag_apart[a][bb].size(); ++bbb){
                                                if(locus_tag_apart[a][aa][aaa] == locus_tag_apart[a][bb][bbb] && prod_acc_flag == 0){ // only prod_acc can be duplicated within and between genomes, locus tag is unique across genomes
                                                    if(aa!=bb){
                                                        if(diag > 1)std::cout << "INSERTING" << '\t' << locus_tag_apart[a][bb][bbb] << '\t' << "onto" << '\t' << x << '\t' << "at pos" << '\t' << pos << '\n';
                                                        final_answer[x][pos].push_back(locus_tag_apart[a][bb][bbb]);
                                                        do_assign = 0;
                                                        goto try_assign;
                                                    }
                                                }
                                                for(int yy = 0; yy < final_answer[x].size(); ++yy){
                                                    for(int yyy = 0; yyy < final_answer[x][yy].size(); ++yyy){
                                                        if(locus_tag_apart[a][bb][bbb] == final_answer[x][yy][yyy]){ // ortholog pair is already contained in final answer, do not duplicate
                                                            do_assign = 0;
                                                            goto next_bbb;
                                                        }else if(diag > 1){
                                                            std::cout << "NO MATCH INNER" << '\t' << locus_tag_apart[a][bb][bbb] << '\t' << final_answer[x][yy][yyy] << '\n';
                                                        }
                                                    }
                                                }
                                                // do insert
                                                if(diag > 1)std::cout << "INSERTING" << '\t' << locus_tag_apart[a][bb][bbb] << '\t' << "onto" << '\t' << x << '\t' << "at pos" << '\t' << pos << '\n';
                                                final_answer[x][pos].push_back(locus_tag_apart[a][bb][bbb]);
                                                do_assign = 0;
                                                next_bbb:;
                                            }
                                        }
                                    }else if(diag > 1){
                                        std::cout << "NO MATCH" << '\t' << locus_tag_apart[a][aa][aaa] << '\t' << final_answer[x][xx][xxx] << '\n';
                                    }
                                }
                            }
                        }
                    }
                }
                if(do_assign == 1){ // If A-B is already in final answer but A-C does not exist, check if B-C exists and if so, insert C onto A-B
                    if(diag>1)std::cout<< "STARTING LOCUS TAG AND FINAL ANSWER COMPARE" << '\n';
                    for(int aa = 0; aa < locus_tag_apart[a].size(); ++aa){
                        int pos = locus_ft[a][aa];
                        for(int aaa = 0; aaa < locus_tag_apart[a][aa].size(); ++aaa){
                            for(int b = 0; b < locus_tag_apart.size(); ++b){
                                if(a == b)continue;
                                for(int bb = 0; bb < locus_tag_apart[b].size();++bb){
                                    for(int bbb = 0; bbb < locus_tag_apart[b][bb].size(); ++bbb){
                                        if(locus_tag_apart[a][aa][aaa] == locus_tag_apart[b][bb][bbb]){
                                            if(diag>1)std::cout << "MATCH BETWEEN LOCUS TAG APART AND LOCUS TAG APART" << '\t' << locus_tag_apart[a][aa][aaa] << '\n';
                                            int cc{0};
                                            if(bb == 0)cc = 1;
                                            for(int ccc = 0; ccc < locus_tag_apart[b][cc].size(); ++ccc){
                                                for(int x = 0; x < final_answer.size(); ++x){
                                                    for(int xx = 0; xx < final_answer[x].size(); ++xx){
                                                        for(int xxx = 0; xxx < final_answer[x][xx].size(); ++xxx){
                                                            if(locus_tag_apart[b][cc][ccc] == final_answer[x][xx][xxx]){
                                                                if(diag>1)std::cout<< "MATCH BETWEEN LOCUS TAG APART AND FINAL ANSWER" << '\t' << locus_tag_apart[b][cc][ccc] << '\n';
                                                                if(diag>1)std::cout << "INSERTING" << '\t' << locus_tag_apart[a][aa][aaa] << '\t' << "onto" << '\t' << x << '\t' << "at pos" << '\t' << pos << '\n';
                                                                final_answer[x][pos].push_back(locus_tag_apart[a][aa][aaa]);
                                                                do_assign = 0;
                                                                goto try_assign;
                                                            }else{
                                                                if(diag > 1)std::cout << "NO MATCH" << '\t' << locus_tag_apart[b][cc][ccc] << '\t' << final_answer[x][xx][xxx] << '\n';
                                                                for(int d = 0; d < locus_tag_apart.size(); ++d){
                                                                    if(b == d)continue;
                                                                    for(int dd = 0; dd < locus_tag_apart[d].size(); ++dd){
                                                                        for(int ddd = 0; ddd < locus_tag_apart[d][dd].size(); ++ddd){
                                                                            if(locus_tag_apart[b][cc][ccc] == locus_tag_apart[d][dd][ddd]){
                                                                                int ee{0};
                                                                                if(d == 0)ee = 1;
                                                                                for(int eee = 0; eee < locus_tag_apart[d][ee].size(); ++eee){
                                                                                    if(locus_tag_apart[d][ee][eee] == final_answer[x][xx][xxx]){
                                                                                        if(diag>1)std::cout<< "MATCH BETWEEN LOCUS TAG APART AND LOCUS TAG APART" << '\t' << locus_tag_apart[b][cc][ccc] << '\t' << "AND MATCH BETWEEN LOCUS TAG APART AND FINAL ANSWER" << '\t' << locus_tag_apart[d][ee][eee] << '\n';
                                                                                        if(diag>1)std::cout << "INSERTING" << '\t' << locus_tag_apart[a][aa][aaa] << '\t' << "onto" << '\t' << x << '\t' << "at pos" << '\t' << pos << '\n';
                                                                                        final_answer[x][pos].push_back(locus_tag_apart[a][aa][aaa]);
                                                                                        do_assign = 0;
                                                                                        goto try_assign;
                                                                                    }
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                try_assign:;
                if(do_assign == 1){ // If A-C could not match A or C to anything in final answer, assign A-C into blank spot on final answer
                    int save_x{-1};
                    for(int x = 0; x < final_answer.size(); ++x){
                        for(int xx = 0; xx < final_answer[x].size(); ++xx){
                            for(int xxx = 0; xxx < final_answer[x][xx].size(); ++xxx){
                                if(final_answer[x][xx][xxx] != ""){
                                    goto next_x;
                                }
                            }
                            
                            
                        }
                        save_x = x;
                        goto assign;
                        next_x:;
                    }
                    assign:;
                    for(int bb = 0; bb < locus_tag_apart[a].size(); ++bb){
                        for(int bbb = 0; bbb < locus_tag_apart[a][bb].size(); ++bbb){
                            int pos = locus_ft[a][bb];
                            if(diag > 1)std::cout << "ASSIGN" << '\t' << locus_tag_apart[a][bb][bbb] << '\t' << "onto" << '\t' << save_x << '\t' << "at pos" << pos << '\n';
                            final_answer[save_x][pos].push_back(locus_tag_apart[a][bb][bbb]);
                        }
                    }
                }
            }

            if(diag > 1){
                std::cout << "FINAL SOG GROUPS" << '\n';
                for(int z = 0; z < final_answer.size(); ++z){
                    for(int zz = 0; zz < final_answer[z].size(); ++zz){
                        for(int zzz = 0; zzz < final_answer[z][zz].size(); ++zzz){
                            std::cout << final_answer[z][zz][zzz] << '\t';
                        }
                        std::cout << '\n';
                    }
                    std::cout << '\n';
                }
            }

            for(int aaa = 0; aaa < final_answer.size(); ++aaa){ // inserts comma between paralogs to match OrthoFinder's output
                for(int bbb = 0; bbb < final_answer[aaa].capacity(); ++bbb){
                    if(final_answer[aaa][bbb].size() > 1){
                    for(int ccc = 0; ccc < final_answer[aaa][bbb].size() - 1; ++ccc){
                        final_answer[aaa][bbb][ccc].insert(final_answer[aaa][bbb][ccc].end(), ',');
                        final_answer[aaa][bbb][ccc].insert(final_answer[aaa][bbb][ccc].end(), ' ');
                    }
                    }
                    if(final_answer[aaa][bbb].size() > 0){ // Reinserts the tab
                        final_answer[aaa][bbb][final_answer[aaa][bbb].size() - 1] += '\t';
                    }
                }
            }

            auto stop = std::chrono::high_resolution_clock::now();
            auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
            if(diag > 0){
                std::cout << "TIME BY MATCH" << '\t' << duration.count() << '\t' << "microseconds" << '\n';
            }

            /* Going out of bounds, prob need to match size and capacity to other loops
            if(diag > 1){
                std::cout << "PRINT NEW" << '\t' << locus_tag.size()<< '\n';
                for(int xx = 0; xx < final_answer.size(); ++xx){
                    for(int yy = 0; yy < final_answer[xx].capacity(); ++yy){
                        for(int zz = 0; zz < final_answer[xx][yy].capacity(); ++zz){
                            std::cout << final_answer[xx][yy][zz] << '\t' << "hey" << '\t';
                        }
                    }
                    std::cout << '\n' <<  "break" << '\n';
                }
                std::cout << "222" << '\n';
                for(int u = 0; u < HOG_master[0][j][0].size(); ++u){
                    std::cout << HOG_master[0][j][0][u] << '\t';
                }            
                std::cout << '\n';
            }
            */
            
            bool new_results{0};
            for(int xx = 0; xx < final_answer.size(); ++xx){ // If not printing all orthofinder, check if there was a change from orthofinder's output. If so, print
                int result_count{0};
                for(int yy = 0; yy < final_answer[xx].capacity(); ++yy){
                    for(int zz = 0; zz < final_answer[xx][yy].capacity(); ++zz){
                        if(final_answer[xx][yy][zz] != ""){
                            ++result_count;
                        }
                    }
                }
                if(result_count == sum && print_all_orthofinder == 0){
                    new_results = 1;
                }else if(result_count != sum && result_count > 0 && print_all_orthofinder == 0){
                    new_results = 0;
                    goto new_loop;
                }
            }
            new_loop:;

            if(new_results == 0){
                for(int xx = 0; xx < final_answer.size(); ++xx){
                    std::string temp_hog_match = "";
                    for(int yy = 0; yy < final_answer[xx].capacity(); ++yy){
                        if(final_answer[xx][yy].size() == 0)temp_hog_match += '\t'; // insert tabs for genomes not contributing a gene to a SOG
                        for(int zz = 0; zz < final_answer[xx][yy].size(); ++zz){
                            if(diag > 1){ 
                                std::cout << "Final answer before inserting to final print array" << '\t' << j << '\t' << xx << '\t' << yy << '\t' << zz << '\t' << final_answer[xx][yy][zz] << '\n';
                            }
                            temp_hog_match += final_answer[xx][yy][zz];
                        }
                    }
                    //std::cout << "temp_hog_match" << '\t' << temp_hog_match << '\n';
                    bool check_only_tab{1};
                    for(int i = 0; i < temp_hog_match.size(); ++i){ // check to make sure the print statement contains more than tabs
                        if(temp_hog_match[i] != '\t'){
                            check_only_tab = 0;
                            break;
                        }
                    }
                    if(check_only_tab == 1){
                        continue;
                    }
                    if(hog_match[j][0].empty() && temp_hog_match != ""){
                        hog_match[j][0] = temp_hog_match;
                    }else if(temp_hog_match != ""){
                        hog_match[j].push_back(temp_hog_match);
                    }
                    if(diag > 0){
                        std::cout << "FINAL SIZE" << '\t' << final_answer.size() << '\n';
                    }
                }   
            }else if(diag > 0){
                std::cout << "New result is false" << '\t' << "HOG" << '\t' << j << '\t' << "size" << '\t' << "size" << '\t' << new_results << '\n';
            }
        }
    }
    if(benchmark == 0){ // removes empty benchmark file as benchmarking was not turned on
        benchmark_file.close();
        char *cstring = new char[out_file_name.length() + 1];
        std::strcpy(cstring, out_file_name.c_str());
        std::remove(cstring);
    }else{
        benchmark_file.close();
    }

    return hog_match;
}

void while_loop(std::vector<int>*** HOG_master, int j, int p, int **zstore, int *z_neg_count, int *z_scaff_count, std::vector<strct_feature_tables_info> *feature_tables_info, int ftk, int m, int window_size, int diag, int ftnk, int pairwise, int combo[][2], int &z_max_operon_count){

    int start_0_counter = 0;
    int samp_back{0};
    int remember_operon{-1};
    int adjust_pre{0};
    int pre_operon_count{0};
    while(start_0_counter != (window_size / 2)){  // loop through genes before gene of interest to calculate the window range
        int n = m - start_0_counter - 1 - samp_back; // -1 so we don't start on m for first loop
        if(n < 0){
            n = (*feature_tables_info)[ftk].HOG.size() + n;  // n is neg so adjust back to positive, wrapping around the genome
        }
        int z = (*feature_tables_info)[ftk].HOG[n];  
        if(z == -999){
            ++samp_back; // go back one more in genome's window as this gene was not assigned to a HOG
            continue;
        }else{
            int y = -1111;       
            for(int i = 0; i < pairwise; ++i){
                if(combo[i][0] == ftk && combo[i][1] == ftnk){ // checks if in pairwise, there is no ortholog between 2 samples when more than 2 samples are present. 
                    y = HOG_master[i][z][1][0];                // may occur if A-B and A-C but not B-C. 
                    break;
                }else if(combo[i][0] == ftnk && combo[i][1] == ftk){
                    y = HOG_master[i][z][0][0];
                    break;
                }
            }
            if((*feature_tables_info)[ftk].scaffold[m] != (*feature_tables_info)[ftk].scaffold[n]){ // if genes are located on different scaffolds, don't cross scaffold
                adjust_pre = ((window_size / 2) - start_0_counter);
                break;
            }else if(y == -9999){
                ++samp_back;
                if(diag > 2){
                    std::cout << "top" << '\t' << "y is -999" << '\t' << y << '\t' << "z is" << '\t' << z << '\t' << "n is" << '\t' << n << '\t' << "m is" << '\t' << m << '\t' << "m at hog ftk" << '\t' << (*feature_tables_info)[ftk].HOG[m] << '\t' << "ftk is" << '\t' << ftk << '\t' << "ftnk is" << '\t' << ftnk << '\t' << "prod acc" << '\t' << (*feature_tables_info)[ftk].prod_acc[n] <<'\n';
                }
            }else if((*feature_tables_info)[ftk].domain != 'e' && (*feature_tables_info)[ftk].operon[n] == remember_operon && (*feature_tables_info)[ftk].operon[n] != -1 && (*feature_tables_info)[ftk].operon[n] != (*feature_tables_info)[ftk].operon[m]){
                ++pre_operon_count; //could track a max within this loop but this will work and take a few extra memory for the wasted space
                ++samp_back;
                if(diag > 2){
                    std::cout << "SAME OPERON PRE" << '\t' << "n" << '\t' << n << '\t' << (*feature_tables_info)[ftk].operon[n] << '\t' << remember_operon << '\t' << ftk << '\t' << "m" << '\t' << m << '\t' << (*feature_tables_info)[ftk].operon[m] << '\t' << (*feature_tables_info)[ftk].locus_tag[n] << '\t' << "samp_back" << '\t' << samp_back <<'\n';
                }
            //}else if(z == (*feature_tables_info)[ftk].HOG[m]){ // same hog as gene of interest, cant double count it
                //++samp_back;
            }else{
                ++start_0_counter;  // there is an ortholog, though not checking if in same window, in both samples
                remember_operon = (*feature_tables_info)[ftk].operon[n];
            }
        }
    }

    start_0_counter = 0;
    int samp_forw{0};  // forward
    remember_operon = -1;
    int adjust_post{0};
    int post_operon_count{0};
    while(start_0_counter != (window_size / 2)){  // same as above but loop through genes after gene of interest to calculate window range
        int n = m + start_0_counter + 1 + samp_forw; // +1 so we don't start on m for first loop
        if(n > ((*feature_tables_info)[ftk].HOG.size() - 1)){
            n = n - (*feature_tables_info)[ftk].HOG.size();  
        }
        int z = (*feature_tables_info)[ftk].HOG[n];  
        if(z == -999){
            ++samp_forw;
            continue;
        }else{
            int y = -1111;        
            for(int i = 0; i < pairwise; ++i){
                if(combo[i][0] == ftk && combo[i][1] == ftnk){
                    y = HOG_master[i][z][1][0];
                    break;
                }else if(combo[i][0] == ftnk && combo[i][1] == ftk){
                    y = HOG_master[i][z][0][0];
                    break;
                }
            }
            if((*feature_tables_info)[ftk].scaffold[m] != (*feature_tables_info)[ftk].scaffold[n]){
                adjust_post = ((window_size / 2) - start_0_counter);
                break;
            }else if(y == -9999){
                ++samp_forw;
                if(diag > 2){
                    std::cout << "bot" << '\t' << "y is -999" << '\t' << y << '\t' << "z is" << '\t' << z << '\t' << "n is" << '\t' << n << '\t' << "m is" << '\t' << m << '\t' << "m at hog ftk" <<    '\t' << (*feature_tables_info)[ftk].HOG[m] << '\t' << "ftk is" << '\t' << ftk << '\t' << "ftnk is" << '\t' << ftnk << '\t' << "prod acc" << '\t' << (*feature_tables_info)[ftk].prod_acc[n] <<'\n';
                }
            }else if((*feature_tables_info)[ftk].domain != 'e' && (*feature_tables_info)[ftk].operon[n] == remember_operon && (*feature_tables_info)[ftk].operon[n] != -1 && (*feature_tables_info)[ftk].operon[n] != (*feature_tables_info)[ftk].operon[m]){
                ++post_operon_count;
                ++samp_forw;
                if(diag > 2){
                    std::cout << "SAME OPERON POST" << '\t' << "n" << '\t' << n << '\t' << (*feature_tables_info)[ftk].operon[n] << '\t' << remember_operon << '\t' << ftk << '\t' << "m" << '\t' << m << '\t' << (*feature_tables_info)[ftk].operon[m] << '\t' << (*feature_tables_info)[ftk].locus_tag[n] << '\t' << "samp_forw" << '\t' << samp_forw <<'\n';
                }
            //}else if(z == (*feature_tables_info)[ftk].HOG[m]){ // same hog as gene of interest, cannot double count
                //++samp_forw;
            }else{
                ++start_0_counter;
                remember_operon = (*feature_tables_info)[ftk].operon[n];
            }
        }
    }
    if(diag > 1){
        std::cout << "prechecking" << '\t' << samp_back << '\t' << samp_forw << std::endl;
        std::cout << "HOG" << '\t' << j << '\t' << "Gene" << '\t' << (*feature_tables_info)[ftk].prod_acc[m] << '\t' << (*feature_tables_info)[ftk].locus_tag[m] << '\t' << "of feature table" << '\t' << (*feature_tables_info)[ftk].name << '\t' << "line #" << '\t' << m << '\n';
        std::cout << "FT_1_line_#" << '\t' << "FT_2_line_#" << '\t' << "HOG_FT_1" << '\t' << "HOG_FT_2" << '\n';
    }

    z_max_operon_count = pre_operon_count + post_operon_count;  
    if(z_max_operon_count == 0){
        z_max_operon_count = 1;
    }

    for(int i = 0; i < window_size + 1; ++i){
        zstore[i] = new int[z_max_operon_count];
        for(int j = 0; j < z_max_operon_count; ++j){
            zstore[i][j] = -999;
        }
    }

    int n{-9999};
    int end{-9999};
    bool adjust_n{0};
    int loop_count = -(window_size / 2);
    int loop_count2{0}; // count from 0 to window_size
    // check if the start of window needs to wrap around geneome or be set to 0
    if(m < ((window_size / 2) + samp_back) - adjust_pre){
        if((*feature_tables_info)[ftk].chromosome == 'c' || (*feature_tables_info)[ftk].chromosome == 'C'){
            int temp = (*feature_tables_info)[ftk].HOG.size();
            n = temp - abs((m - ((window_size / 2) - adjust_pre)) - samp_back);
            adjust_n = 1;
            if(diag > 2){
                std::cout << "UnderFLOW MATH" << '\t' << "m" << '\t' << m <<  '\t' << samp_back << std::endl;
                std::cout << "temp" << '\t' << temp << std::endl;
                std::cout << "n in underflow" << '\t' << n << std::endl;
            }
        }else{ // chromosome is linear
            if(diag > 2){
                std::cout << "Linear n" << '\t' << (*feature_tables_info)[ftk].HOG[m] << '\t' << "ftk" << '\t' << ftk << '\t' << "m" << '\t' << m << '\n';
            }
            n = 0;
        }
    }else{
        n = m - (((window_size / 2) - adjust_pre) + samp_back);
    }
    // check if the end of window needs to wrap around genome or be set to max
    if(m >= (*feature_tables_info)[ftk].HOG.size() - (window_size / 2) - samp_forw){ // might have to include adjust_post like we did adjust_pre on the m < check in the above if
        if(diag > 2){
            std::cout << "OVERFLOW MATH" << '\t' << (*feature_tables_info)[ftk].HOG.size() << '\t' << m << '\t' << samp_forw << '\n';
        }
        if((*feature_tables_info)[ftk].chromosome == 'c' || (*feature_tables_info)[ftk].chromosome == 'C'){

            int temp = (*feature_tables_info)[ftk].HOG.size() - ((m + samp_forw) + ((window_size / 2) - adjust_post));
            end = abs(temp);
            if(adjust_n == 1){
                std::cout << "Error, too large of window size for genome" << '\t' <<  (*feature_tables_info)[ftk].name << '\t' << "Underflowing and OverFlowing the feature table file length." << '\t' << (*feature_tables_info)[ftk].HOG.size()   << '\n';
                std::exit(EXIT_FAILURE);
            }
            adjust_n = 1;
        }else{  // chromosome is linear
            if(diag > 2){
                std::cout << "Linear END" << '\t' << (*feature_tables_info)[ftk].HOG[m] << '\t' << "ftk" << '\t' << ftk << '\t' << "m" << '\t' << m << '\t' << (*feature_tables_info)[ftk].chromosome << '\n';
            }
            end = (*feature_tables_info)[ftk].HOG.size();
        }
    }else{
        end = m + samp_forw + ((window_size / 2) - adjust_post) + 1;                        
    }

    if(diag > 0){
        std::cout << "END" << '\t' << end << '\t' << "n" << '\t' << n << '\t' << samp_forw << '\t' << samp_back << '\t' << (*feature_tables_info)[ftk].chromosome << '\n';
    }

    z_scaff_count[p] = adjust_pre + adjust_post;

    remember_operon = -1;
    int operon_2D{0};
    while(n != end || adjust_n != 0){ // loop through range of window storing the HOG of each gene within window
        int z = (*feature_tables_info)[ftk].HOG[n];
        if(n == m){
            if(diag > 1){
                std::cout << "n in window is" << '\t' << n  << '\t' << z << '\n';
            } 
            ++n;
            if(adjust_n == 1){ // check if n should be set to 0 after completing the wrap around
                if(n >= (*feature_tables_info)[ftk].HOG.size()){
                    if(diag > 2){
                        std::cout << "adjusting n to 0" << '\n';
                    }
                    n = 0;
                    adjust_n = 0;
                }else{
                    if(diag > 2){
                        std::cout << "n not greater than (*feature_tables_info)[ftk].HOG.size()" << '\t' << n << '\t' << (*feature_tables_info)[ftk].HOG.size() << std::endl;
                        std::cout << "n not less than end" << '\t' << n << '\t' << end << std::endl;
                    }   
                }
            }
            continue;
        }
        if(diag > 3){
            std::cout << "n" << '\t' << n << '\t' << "end" << '\t' << end << std::endl;
        }

        int y = -1;
        if(z == -999){ // gene not assigned to HOG
            z_neg_count[p] += 1;
        }else{
            for(int i = 0; i < pairwise; ++i){
                if(combo[i][0] == ftk && combo[i][1] == ftnk){
                    y = HOG_master[i][z][1][0];
                    break;
                }else if(combo[i][0] == ftnk && combo[i][1] == ftk){
                    y = HOG_master[i][z][0][0];
                    break;
                }
            }
            if(y == -9999){

            }else if(y == -1){ // safety check
                std::cout << "y is -1, check code" << '\n';  // prob on ++loopcount2 below corrupting combo[][]
            }else if((*feature_tables_info)[ftk].scaffold[m] != (*feature_tables_info)[ftk].scaffold[n] && z != -999){  // genes on different scaffolds
                z_scaff_count[p] += 1;
                if(diag > 1){
                    std::cout << "Scaff count" << '\t' << z << '\t' << z_scaff_count[p] << '\n'; 
                }
                if(diag > 1){
                    std::cout << "Continue m & n scaffold mismatch" << '\t' << m << '\t' << (*feature_tables_info)[ftk].scaffold[m] << '\t' << n << '\t' << (*feature_tables_info)[ftk].scaffold[n] <<  '\n';
                }
            }else if(remember_operon == (*feature_tables_info)[ftk].operon[n] && remember_operon != -1 && (*feature_tables_info)[ftk].operon[n] != (*feature_tables_info)[ftk].operon[m]){ // if genes are in same operon and we are checking for that
                if(diag > 2){
                    std::cout << n << '\t' << z << '\t' << "BOT SAME OPERON" << '\t' << remember_operon << '\t' << loop_count2 - 1 << '\t' << operon_2D << '\n';
                }
                if(operon_2D > z_max_operon_count){
                    std::cout << "operon_2d greater than z_max_operon_count" << '\n';
                    std::exit(EXIT_FAILURE);
                }else{
                    //std::cout << "2D_operon" << '\t' << operon_2D << '\t' << "max" << '\t' << z_max_operon_count << '\n';
                }
                zstore[loop_count2 - 1][operon_2D] = z; // -1 because the loop_count2 is incremented after storing into zstore in else below. Need to store these operons members into same 1stD array index. so -1
                ++operon_2D;
                
                if(loop_count2 > window_size){
                    std::cout << "loop_count2 bigger than windowsize + 1" << '\n';
                    std::exit(EXIT_FAILURE);
                }else{
                    //std::cout << "loop_count2" << '\t' << loop_count2 << '\n';
                }
            //}else if(z == (*feature_tables_info)[ftk].HOG[m]){ // same hog as gene of interest, cannot double count
                
            }else{ // store the HOG number of genes within window
                remember_operon = (*feature_tables_info)[ftk].operon[n];
                if(loop_count2 > window_size){
                    std::cout << "loop_count2 bigger than windowsize + 1" << '\n';
                    std::exit(EXIT_FAILURE);
                }else{
                    //std::cout << "loop_count2" << '\t' << loop_count2 << '\n';
                }
                zstore[loop_count2][0] = z;
                ++loop_count2;
                operon_2D = 0;
                if(diag > 1){
                    std::cout << n  << '\t' << z << '\t' << (*feature_tables_info)[ftk].operon[n] << '\n';    
                }
                if(diag > 3){
                    std::cout << "Matching scaffolds" << '\t'<< (*feature_tables_info)[ftk].scaffold[m] << '\t' << (*feature_tables_info)[ftk].scaffold[n] << '\n';
                }    
            }
        }
        
        ++n;
        if(adjust_n == 1){
            if(n >= (*feature_tables_info)[ftk].HOG.size()){
                if(diag > 2){
                    std::cout << "adjusting n to 0" << '\n';
                }
                n = 0;
                adjust_n = 0;
            }else{
                if(diag > 2){
                    std::cout << "n not greater than (*feature_tables_info)[ftk].HOG.size()" << '\t' << n << '\t' << (*feature_tables_info)[ftk].HOG.size() << std::endl;
                    std::cout << "n not less than end" << '\t' << n << '\t' << end << std::endl;
                }
            }
        }
    }
}

void find_substr(bool &found_t, size_t &pos, std::vector<std::vector<std::string>> locus_tag, int xx, int yy, std::string &search){

    std::string temp = locus_tag[xx][yy]; // locus_tag will also be final_answer
    if(temp.find('\t', pos + 1) != std::string::npos){ // if tab found, store from pos to final char before tab
        if(temp.find('\t', temp.find('\t', pos + 1)) && found_t == 1){ // if a tab was found and more tabs will be found
            search = temp.substr(pos + 1, temp.find('\t'));
        }else{ // no more tabs after this tab
            search = temp.substr(pos, temp.find('\t'));
        }
        found_t = 1;
    }else if(found_t == 1){
        search = temp.substr(pos + 1, std::string::npos); //store from last tab to end of string
    }else{
        search = temp.substr(pos, std::string::npos);  // no more tab, store rest of string. If not tab ever, store entire string
    }
    pos = temp.find('\t', pos + 1); // set pos_a to next tab or end of string
}
