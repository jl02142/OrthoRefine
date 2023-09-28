#include <iostream>
#include <fstream>
#include <sstream>


int main(int argc, char* argv[]){

    std::cout << argv[1] << '\t' << argv[2] << '\n';

    std::string infile_name_prot = argv[2];
    std::ifstream infile_prot(infile_name_prot);
    if(infile_prot){
        
    }else{
        std::cout << "Error: file " << infile_name_prot << " not found" << "\n";
        std::exit(EXIT_FAILURE);
    }

    std::string store_string;
    std::string temp_string = "";    
    bool h_flag{0}; // header flag
    std::stringstream ss;
    char c;


    store_string = "";
    infile_prot.seekg(0, std::ios::end);     // Go to end of file
    store_string.resize(infile_prot.tellg()); // resize store_string to size of infile
    infile_prot.seekg(0, std::ios::beg);     // Go to start of file
    infile_prot.read(&store_string[0], store_string.size());  // Read file into store_string
    infile_prot.close();

    temp_string = "";
    h_flag = 0;
    int prot_file_line_count{0};


    ss.clear();
    ss << store_string;
    while(ss >> std::noskipws >> c){
        if(h_flag == 1 && c == '\n'){
            ++prot_file_line_count;
            temp_string = "";
            h_flag = 0;
        }else if(h_flag == 0 && temp_string == "CDS"){
            h_flag = 1;
        }else if(h_flag == 0 && c == ' '){
            temp_string = "";
        }else if(h_flag == 0 && c == '\n'){
            temp_string = "";
        }else if(h_flag == 0){
            temp_string += c;
        }
    }

    std::cout << "prot_line count" << '\t' << prot_file_line_count << '\n';
    // prot_file_line_count should be equal to cat $file | grep -c " CDS "

    std::string **original_prot = new std::string *[prot_file_line_count];
    for(int i = 0; i < prot_file_line_count; ++i){
        original_prot[i] = new std::string[2];
        for(int j = 0; j < 2; ++j){
            original_prot[i][j] = "";
        }
    }

    ss.clear();
    ss << store_string;
    h_flag = 0;
    bool GN_flag{0};
    bool prestore_uniprot{0};
    bool uniprot_flag{0};
    int temp_prot_line_count{0};
    while(ss >> std::noskipws >> c){
        if(c == '\n'){
            temp_string = "";
        }else if(GN_flag == 1 && c == '"'){
            GN_flag = 0;
            int i = temp_prot_line_count;
            
            original_prot[i][0] = temp_string;
            temp_string = "";
        }else if(uniprot_flag == 1 && c == '"'){
            uniprot_flag = 0;
            h_flag = 0;
            int i = temp_prot_line_count;
            original_prot[i][1] = temp_string;
            ++temp_prot_line_count;
            temp_string = "";
        }else if(h_flag == 1 && temp_string == "/locus_tag=\""){
            temp_string = "";
            temp_string += c;
            GN_flag = 1;
        }else if(prestore_uniprot == 1 && c == ':'){
            temp_string = "";
            uniprot_flag = 1;
            prestore_uniprot = 0;            
        }else if(h_flag == 1 && temp_string == "/db_xref=\"UniProtKB/"){
            prestore_uniprot = 1;
        }else if(h_flag == 0 && temp_string == "CDS"){
            h_flag = 1;
        }else if(c == ' '){
            temp_string = "";
        }else{
            temp_string += c;
        }
    }

    /*
    for(int i = 0; i < prot_file_line_count; ++i){
        std::cout << i << '\t';
        for(int j = 0; j < 2; ++j){
            std::cout << original_prot[i][j] << '\t';
        }
        std::cout << '\n';
    }
    */


    std::string infile_name_NCBI = argv[1];
    std::ifstream infile_NCBI(infile_name_NCBI);
    if(infile_NCBI){
        
    }else{
        std::cout << "Error: file " << infile_name_NCBI << " not found" << "\n";
        std::exit(EXIT_FAILURE);
    }

    h_flag = 0;
    ss.clear();
    store_string = "";
    infile_NCBI.seekg(0, std::ios::end);     // Go to end of file
    store_string.resize(infile_NCBI.tellg()); // resize store_string to size of infile
    infile_NCBI.seekg(0, std::ios::beg);     // Go to start of file
    infile_NCBI.read(&store_string[0], store_string.size());  // Read file into store_string
    infile_NCBI.close();

    int NCBI_file_line_count{0};
    int NCBI_file_tab_count{1};

    ss << store_string;
    while(ss >> std::noskipws >> c){
        if(c == '\n'){
            h_flag = 1;
            ++NCBI_file_line_count;
        }else if(h_flag == 0 && c == '\t'){
            ++NCBI_file_tab_count;
        }
    }

    std::cout << "NCBI_line count" << '\t' << NCBI_file_line_count << '\n';
    std::cout << "NCBI_tab_count" << '\t' << NCBI_file_tab_count << '\n';

    std::string **original_NCBI = new std::string*[NCBI_file_line_count];
    for(int i = 0; i < NCBI_file_line_count; ++i){
        original_NCBI[i] = new std::string[NCBI_file_tab_count];
        for(int j = 0; j < NCBI_file_tab_count; ++j){
            original_NCBI[i][j] = "";
        }
    }

    int temp_tab_count{0};
    int temp_line_count{0};
    ss.clear();
    ss << store_string;
    
    while(ss >> std::noskipws >> c){
        if(c == '\n'){
            int i = temp_line_count;
            int j = temp_tab_count;
            original_NCBI[i][j] = temp_string;
            temp_string = "";
            temp_tab_count = 0;
            ++temp_line_count;
        }else if(c == '\t'){
            int i = temp_line_count;
            int j = temp_tab_count;
            original_NCBI[i][j] = temp_string;
            temp_string = "";
            ++temp_tab_count;
        }else if(c != '\t'){
            temp_string += c;
        }
    }

    /*
    for(int i = 0; i < NCBI_file_line_count; ++i){
        std::cout << i << '\t';
        for(int j = 0; j < NCBI_file_tab_count; ++j){
            std::cout << original_NCBI[i][j] << '\t';
        }
        std::cout << '\n';
    }
    */

    // loops below are inefficent when there are multiple old locus tags to check but o well
    for(int i = 0; i < prot_file_line_count; ++i){
        for(int j = 0; j < NCBI_file_line_count; ++j){
            if(original_NCBI[j][0] != "CDS"){
                continue;
            }
            if(original_prot[i][0] == original_NCBI[j][16]){
                original_NCBI[j][16] = original_prot[i][1];
                goto next_i;
            }
        }
        for(int j = 0; j < NCBI_file_line_count; ++j){
            if(original_NCBI[j][0] != "gene"){  
                continue;
            }
            size_t found = original_NCBI[j][19].find("old_locus_tag=");
            if(found != std::string::npos){
                std::string old_locus_tag = original_NCBI[j][19].substr(found + 14);
                if(original_prot[i][0] == old_locus_tag){
                    found = original_NCBI[j][16].find(original_NCBI[j + 1][16]);
                    if(found != std::string::npos){
                        original_NCBI[j + 1][16] = original_prot[i][1];
                    }else{
                        std::cout << "Warning: " << '\t' << old_locus_tag << " found on gene line but does not match following CDS line" << '\n';
                    }
                    goto next_i;
                }
                size_t pos = 0;
                size_t found_comma;
                do{
                    //std::cout << "old locus tag" << '\t' << old_locus_tag << '\t' << pos << '\n';
                    found_comma =  old_locus_tag.substr(pos).find(",");
                    if(found_comma == std::string::npos){
                        std::string compare = old_locus_tag.substr(pos);
                        //std::cout << "compare top" << '\t' << compare << '\n';
                        if(original_prot[i][0] == compare){
                            found = original_NCBI[j][16].find(original_NCBI[j + 1][16]);
                            if(found != std::string::npos){
                                original_NCBI[j + 1][16] = original_prot[i][1];
                            }else{
                                std::cout << "Warning: " << '\t' << old_locus_tag << " found on gene line but does not match following CDS line" << '\n';
                            }
                            goto next_i;
                        }
                    }
                    //std::cout << i << '\t' << j << '\t' << found_comma << '\n';
                    if(found_comma != std::string::npos){
                        std::string compare = old_locus_tag.substr(pos, found_comma);
                        //std::cout << "compare bot" << '\t' << compare << '\n';
                        pos = pos + found_comma + 1;
                        if(original_prot[i][0] == compare){
                            found = original_NCBI[j][16].find(original_NCBI[j + 1][16]);
                            if(found != std::string::npos){
                                original_NCBI[j + 1][16] = original_prot[i][1];
                            }else{
                                std::cout << "Warning: " << '\t' << old_locus_tag << " found on gene line but does not match following CDS line" << '\n';
                            }
                            goto next_i;
                        }
                    }
                }while(found_comma != std::string::npos);
            }
        }
        next_i:;
    }

    /*
    for(int i = 0; i < NCBI_file_line_count; ++i){
        size_t found = original_NCBI[i][19].find("old_locus_tag=");
        if(found != std::string::npos){
            std::cout << "found" << '\t' << original_NCBI[i][19] << '\n';
            found = original_NCBI[i][16].find(original_NCBI[i + 1][16]);
            if(found != std::string::npos){
                std::cout << "replace" << '\t' << original_NCBI[i + 1][16] << '\n';
                std::string temp = original_NCBI[i][19].substr(found + 14);
                original_NCBI[i + 1][16] = temp;
                std::cout << "DONE" << '\t' << original_NCBI[i + 1][16] << '\t' << temp << '\n';
            }else{
                //std::cout << "not_found" << '\t' << original_NCBI[i][19] << '\t' << original_NCBI[i + 1][16] << '\n';
            }
        }
    }
    */


    /*
    for(int i = 0; i < NCBI_file_line_count; ++i){
        std::cout << i << '\t';
        for(int j = 0; j < NCBI_file_tab_count; ++j){
            std::cout << original_NCBI[i][j] << '\t';
        }
        std::cout << '\n';
    }
    */


    
    /*    
    for(int i = 0; i < NCBI_file_line_count; ++i){
        if(original_NCBI[i][0] == "CDS"){
            for(int j = 0; j < prot_file_line_count; ++j){
                if(original_NCBI[i][16] == original_prot[j][0]){
                    //std::cout << "MATCH" << '\n';
                    original_NCBI[i][16] = original_prot[j][1];
                }else{
                    //std::cout << "NOT_MATCH" << '\t' << original_NCBI[i][16] << '\t' << original_prot[j][0] << '\n';
                }
            }
        }
    }
    */


    /* 
    for(int i = 0; i < NCBI_file_line_count; ++i){
        //std::cout << i << '\t';
        for(int j = 0; j < NCBI_file_tab_count; ++j){
            std::cout << original_NCBI[i][j] << '\t';
        }
        std::cout << '\n';
    }
    */



    std::string outfile_name = infile_name_NCBI;
    std::size_t found = infile_name_NCBI.find(".");  // find first period
    found = infile_name_NCBI.find(".", found + 1);  // find second period
    std::string replace = "_locus_tag_replaced_with_uniprot.txt";
    outfile_name.replace(infile_name_NCBI.find(".", found), replace.length(), replace);
    
    std::ofstream outfile(outfile_name);
    for(int i = 0; i < NCBI_file_line_count; ++i){
        for(int j = 0; j < NCBI_file_tab_count; ++j){
            outfile << original_NCBI[i][j] << '\t';
        }
        outfile << '\n';
    }




}
