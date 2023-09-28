#include <iostream>
#include <fstream>
#include <sstream>


int main(int argc, char* argv[]){

    int ortho_line_count{0};
    int ortho_tab_count{1};
    int ortho_comma_max_count{0};
    int temp_ortho_comma_count{0};
    bool h_flag{0}; // header flag

    std::string ortho_infile_name = argv[1];
    std::ifstream ortho_infile(ortho_infile_name);
    if(ortho_infile){
        
    }else{
        std::cout << "Error: file " << ortho_infile_name << " not found" << "\n";
        std::exit(EXIT_FAILURE);
    }

    std::string store_string = "";
    std::string temp_string = "";    
    std::stringstream ss;
    char c;

    ortho_infile.seekg(0, std::ios::end);     // Go to end of file
    store_string.resize(ortho_infile.tellg()); // resize store_string to size of infile
    ortho_infile.seekg(0, std::ios::beg);     // Go to start of file
    ortho_infile.read(&store_string[0], store_string.size());  // Read file into store_string
    ortho_infile.close();

    ss << store_string;
    while(ss >> std::noskipws >> c){
        if(h_flag == 0 && c == '\t'){
            ++ortho_tab_count;
        }else if(h_flag == 1 && c == '\n'){
            ++ortho_line_count;
            if(temp_ortho_comma_count > ortho_comma_max_count){
                ortho_comma_max_count = temp_ortho_comma_count;
            }
            temp_ortho_comma_count = 0;
        }else if(h_flag == 1 && c == ','){
            ++temp_ortho_comma_count;
        }else if(h_flag == 0 && c == '\n'){
            h_flag = 1;
        }
    }

    //std::cout << "Ortho line count" << '\t' << ortho_line_count << '\n';
    //std::cout << "Ortho tab count" << '\t' << ortho_tab_count << '\n';

    std::string ***ortho_store = new std::string **[ortho_line_count];
    for(int i = 0; i < ortho_line_count; ++i){
        ortho_store[i] = new std::string*[ortho_tab_count];
        for(int j = 0; j < ortho_tab_count; ++j){
            ortho_store[i][j] = new std::string[ortho_comma_max_count];
            for(int k = 0; k < ortho_comma_max_count; ++k){
                ortho_store[i][j][k] = "";
            }
        }
    }

    ss.clear();
    ss << store_string;
    int temp_ortho_line_count{0};
    int temp_ortho_tab_count{0};
    temp_ortho_comma_count = 0;
    h_flag = 0;
    while(ss >> std::noskipws >> c){
        if(h_flag == 1 && c == '\t'){
            int i = temp_ortho_line_count;
            int j = temp_ortho_tab_count;
            int k = temp_ortho_comma_count;
            ortho_store[i][j][k] = temp_string;
            ++temp_ortho_tab_count;
            temp_ortho_comma_count = 0;
            temp_string = "";
        }else if(h_flag == 1 && c == '\n'){
            int i = temp_ortho_line_count;
            int j = temp_ortho_tab_count;
            int k = temp_ortho_comma_count;
            ortho_store[i][j][k] = temp_string;
            ++temp_ortho_line_count;
            temp_ortho_tab_count = 0;
            temp_ortho_comma_count = 0;
            temp_string = "";
        }else if(h_flag == 1 && c == ','){
            int i = temp_ortho_line_count;
            int j = temp_ortho_tab_count;
            int k = temp_ortho_comma_count;
            ortho_store[i][j][k] = temp_string;
            temp_string = "";
            ++temp_ortho_comma_count;
        }else if(h_flag == 0 && c == '\n'){
            h_flag = 1;
        }else if(c != ' '){
            temp_string += c;
        }
    }

    for(int i = 0; i < ortho_line_count; ++i){
        for(int j = 3; j < ortho_tab_count; ++j){
            for(int k = 0; k < ortho_comma_max_count; ++k){
                if(ortho_store[i][j][k] == ""){
                    continue;
                }
                for(int m = j + 1; m < ortho_tab_count; ++m){
                    for(int n = 0; n < ortho_comma_max_count; ++n){
                        if(ortho_store[i][m][n] == ""){
                            continue;
                        }
                        std::cout << ortho_store[i][j][k] << '\t' << ortho_store[i][m][n] << '\n';
                    }
                }
            }
        }
    }
}
