// script to convert orthoxml (OMA) to orthofinder (tsv) format for use with OrthoRefine
// J Ludwig. 2024.

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <filesystem> // std::filesystem::directory_iterator

int main(int argc, char* argv[]){
    
    std::cout << "This script was designed for use with OrthoRefine. The orthoxml file standard is very inconsistent (geneId vs protId, which one is it?!?). Please check the output file for errors." << '\n';
    std::cout << "Reminder: the input file with the GCF accession must have the same order as the Orthoxml file (samples appear in the same order) if using with OrthoRefine" << '\n'; // OrthoRefine takes order from the input file
    
    std::string inputFileName = "";
    std::string outputFileName = "";
    std::string orthoxmlHogFileName = "";
    std::string orthoxmlGroupFileName = "";
    for(int i = 1; i < argc; ++i){
        std::string arg = argv[i];
        if(arg == "--input"){
            inputFileName = argv[++i];
        }else if (arg == "--output"){
            outputFileName = argv[++i];
        }else if (arg == "--orthoxmlHog"){
            orthoxmlHogFileName = argv[++i];
        }else if(arg == "--orthoxmlGroup"){
            orthoxmlGroupFileName = argv[++i];
        }
        else if(arg == "--help" || arg == "-h"){
            std::cout << "Usage: ./orthoxml_convert.exe --input <input file> --output <output file> --orthoxmlHog <orthoxmlHog file> --orthoGroup <orthoxmlGroup file>" << '\n' << "If input not provdied, search orthoxml file for GCF or GCA accession, download them, and use them as input" << '\n'; 
            std::exit(EXIT_SUCCESS);
        }else{
            std::cout << "Error: incorrect option provided" << '\t' << argv[i] << '\n';
            std::cout << "Usage: ./orthoxml_convert.exe --input <input file> --output <output file> --orthoxml <orthoxmlHog file> --orthoGroup <orthoxmlGroup file>" << '\n' << "If input not provdied, search orthoxml file for GCF or GCA accession, download them, and use them as input" << '\n'; 
            std::exit(EXIT_FAILURE);
        }
    }
    if(outputFileName == "" || orthoxmlHogFileName == "" || orthoxmlGroupFileName == ""){
        std::cout << "Error: missing output or orthoxml file name" << '\n';
        std::exit(EXIT_FAILURE);
    }

    std::string line;
    std::vector<std::string> gcX;  // GCF or GCA accession numbers
    std::ifstream orthoxml;
    orthoxml.open(orthoxmlHogFileName); // orthoxml file
    while(std::getline(orthoxml, line)){
            if(line.find("database name=") != std::string::npos){
                size_t pos;
                if(line.find("GCA_") != std::string::npos){
                    pos = line.find("GCA_"); // set pos to the end of "name="
                }else if(line.find("GCF_") != std::string::npos){
                    pos = line.find("GCF_");
                }else{
                    gcX.push_back("NONE"); // no GCF or GCA found
                    continue;
                }
                size_t pos2 = line.find_first_not_of("0123456789.", pos + 4); // set pos2 to the quotation mark after the GCA or GCF
                std::string temp = line.substr(pos, pos2 - pos); // get the name between pos and pos2
                gcX.push_back(temp);
            }
    }
    orthoxml.close();
    if(inputFileName == ""){
        std::cout << "No input file provided. Will find GCA_ and GCF_ in database." << '\n';
        if(!std::filesystem::exists("download_ft_fafiles.sh")){
            std::cout << "Error: download_ft_fafiles.sh not found in current directory" << '\n';
            std::exit(EXIT_FAILURE);
        }
        std::string inOutputFile = outputFileName + "_orthorefine_input.txt";
        std::ofstream output;
        output.open(inOutputFile);
        for(int i = 0; i < gcX.size(); ++i){
            if(gcX[i] != "NONE"){
                output << gcX[i] << '\n';
            }
        }
        output.close();
        std::string command = "./download_ft_fafiles.sh " + inOutputFile; // use download_ft_fafiles.sh to download the .faa and feature table files
        system(command.c_str());
        inputFileName = outputFileName + "_orthorefine_input.txt";
    }
        
    std::ifstream inputNames;  // input file with GCF accession numbers
    inputNames.open(inputFileName);
    std::vector<std::string> fileNames;
    while(std::getline(inputNames, line)){
        fileNames.push_back(line);
    }
    inputNames.close();
    //for(int i = 0; i < fileNames.size(); ++i){
    //    std::cout << "filenames" << '\t' << fileNames[i] << '\n';
    //}
    std::vector<std::string> fullInputNames;
    std::string path = std::filesystem::current_path();
    for(int i = 0; i < fileNames.size(); ++i){ // search the current directory for the full file name
        for(const auto &file : std::filesystem::directory_iterator(path)){
            if(file.path().extension() != ".faa"){ // file must end with .faa
                continue;
            }
            size_t found = file.path().string().find(fileNames[i]);
            if(found != std::string::npos){
                for(const auto &file2 : std::filesystem::directory_iterator(path)){
                    if(file2.path().string().find("_feature_table.txt")){ 
                        size_t found2 = file2.path().string().find(fileNames[i]);
                        if(found2 != std::string::npos){ // require matching feature table file
                            std::string filename = file.path().filename().string();
                            fullInputNames.push_back(filename.substr(0, filename.length() - 4)); // strip the ".faa" to match OrthoFinder output
                            break;
                        }
                    }
                }
                break;
            }
        }
    }
    //for(int i = 0; i < fullInputNames.size(); ++i){
    //    std::cout << i << '\t' << "fullinputnames" << '\t' << fullInputNames[i] << '\n';
    //}

    std::ifstream orthoxmlGroup;
    orthoxmlGroup.open(orthoxmlGroupFileName); // orthoxml Group file
    orthoxml.open(orthoxmlHogFileName); // orthoxml HOG file
    std::vector<std::vector<std::vector<std::string>>> orthoGroup;
    std::vector<int> gene_id;
    std::vector<std::string> geneId;
    std::vector<int> genomeLastNumb; // which genome the gene belongs to by last gene_ID number
    std::vector<int> genomeFirstNumb; // which genome the gene belongs to by first gene_ID number
    std::vector<std::vector<std::string>> temp_vector;
    temp_vector.resize(fullInputNames.size());
    bool inDB{0}; // checks if sample has GCF or GCA
    std::string line2 = "";
    while(std::getline(orthoxml, line)){
        std::string gcXName = "";
        if(line.find("database name=") != std::string::npos){
            size_t pos;
            if(line.find("GCA_") != std::string::npos){
                pos = line.find("GCA_"); // set pos to the end of "name="
            }else if(line.find("GCF_") != std::string::npos){
                pos = line.find("GCF_");
            }else{
                inDB = 0;
                continue;
            }
            size_t pos2 = line.find_first_not_of("0123456789.", pos + 4); // set pos2 to the character after the GCA or GCF accession number
            gcXName = line.substr(pos, pos2 - pos); // get the name between pos and pos2
            for(int i = 0; i < fullInputNames.size(); ++i){
                std::string fullFileName = fullInputNames[i];
                if(fullFileName.find(gcXName) != std::string::npos){
                    //std::cout << "FOUND" << '\t' << fullInputNames[i] << std::endl;
                    inDB = 1;
                    break;
                }
            }
            if(inDB == 0){
                //std::cout << "NOT FOUND" << '\t' << gcXName << '\n';
            }
        }else if(line.find("<genes>") != std::string::npos && inDB == 1){
            int i = gene_id.size(); // not -1 because minus 1 points to previous <genes> line and we need the current <genes> line we are about to process
            bool eof_flag{0}; // only permits the reset of file to start once per loop
            while(std::getline(orthoxml, line) && line.find("</genes>") == std::string::npos){
                size_t pos = line.find("<gene id=") + 10; // set pos to the end of "<gene id"
                size_t pos2 = line.find('"', pos); // set pos2 to the quotation mark after the gene id
                std::string temp = line.substr(pos, pos2 - pos); // get the gene id between pos and pos2
                gene_id.push_back(std::stoi(temp));
                pos = line.find("geneId=");
                pos2 = line.find('"', pos);
                std::string GCX2 = "";
                if(pos == std::string::npos){
                    pos = line.find("protId="); // some orthoxml files use protId instead of geneId
                    if(pos == std::string::npos){
                        continue;
                    }else{
                        pos += 8; // set pos to the end of "protId="
                    }
                    std::string protTemp = line.substr(pos, pos2 - pos);
                    size_t posGroup;
                    size_t posGroup2;
                    while(std::getline(orthoxmlGroup, line2) && line2.find("<groups>") == std::string::npos){
                        if(line2.find("database name=") != std::string::npos){
                            if(line2.find("GCA_") != std::string::npos){
                                posGroup = line2.find("GCA_"); // set pos to the end of "name="
                                posGroup2 = line2.find_first_not_of("0123456789.", posGroup + 4); // set pos2 to the character after the GCA or GCF accession number
                                GCX2 = line2.substr(posGroup, posGroup2 - posGroup); // get the name between pos and pos2
                            }else if(line2.find("GCF_") != std::string::npos){
                                posGroup = line2.find("GCF_");
                                posGroup2 = line2.find_first_not_of("0123456789.", posGroup + 4); // set pos2 to the character after the GCA or GCF accession number
                                GCX2 = line2.substr(posGroup, posGroup2 - posGroup); // get the name between pos and pos2
                            }else{
                                continue;
                            }
                        }
                        if(GCX2 != gcXName){ //GCX from oma-group must match to GCX from oma-hog
                            continue;
                        }else{
                            do{
                                size_t posProt2 = line2.find("protId=");
                                size_t posProt3 = line2.find('"', posProt2);
                                if(posProt2 == std::string::npos){
                                    continue;
                                }else{
                                    posProt2 += 8; // set pos to the end of "protId="
                                } 
                                if(line2.substr(posProt2, posProt3 - posProt2) == protTemp){
                                    size_t posGeneId = line2.find("geneId=");
                                    if(posGeneId == std::string::npos){
                                        continue;
                                    }else{
                                        posGeneId += 8; // set pos to the end of "geneId="
                                    }
                                    size_t posGeneId2 = line2.find('"', posGeneId);
                                    temp = line2.substr(posGeneId, posGeneId2 - posGeneId);
                                    goto needbreak2;
                                }
                            }while(getline(orthoxmlGroup, line2) && line2.find("</genes>") == std::string::npos);
                        }
                        // if orthoxmlGroup reaches end of file, reset to beginning of file to search again as this search remembers its last location
                        if(orthoxmlGroup.eof() && eof_flag == 0){
                            orthoxmlGroup.clear();
                            orthoxmlGroup.seekg(0, std::ios::beg);
                            eof_flag = 1;
                        }
                    }
                    needbreak2:;
                }else{
                    pos += 8; // set pos to the end of "geneId="
                    temp = line.substr(pos, pos2 - pos);
                }
                geneId.push_back(temp);
            }
            genomeFirstNumb.push_back(gene_id[i]); // first gene id for each genome, used to determine which genome the gene belongs to in below loop
            genomeLastNumb.push_back(gene_id.back());  // last gene id for each genome, used to determine which genome the gene belongs to in below loop
        }else if(line.find("<genes>") != std::string::npos && inDB == 0){ // handle missing data from database
            while(std::getline(orthoxml, line) && line.find("</genes>") == std::string::npos){
                gene_id.push_back(-1);
                geneId.push_back("");
            }
        }else if(line.find("<orthologGroup") != std::string::npos){
            
            while(std::getline(orthoxml, line)){
                //std::cout << line << std::endl;
                if(line.find("<paralogGroup>") != std::string::npos || line.find("</paralogGroup>") != std::string::npos){
                    //std::cout << "PARA" << '\n';
                    continue;
                }else if(line.find("<orthologGroup") != std::string::npos){
                    //std::cout << "ORTHO" << '\n';
                    continue;
                }else if(line.find("</orthologGroup>") != std::string::npos){
                    for(int i = 0; i < temp_vector.size(); ++i){
                        for(int j = 0; j < temp_vector[i].size(); ++j){
                            if(temp_vector[i][j] != ""){
                                //std::cout << "PUSHING ORTHOGROUP at pos" << '\t' << orthoGroup.size() + 1 << '\n';
                                orthoGroup.push_back(temp_vector);
                                temp_vector.clear();
                                temp_vector.resize(fullInputNames.size()); // have to resize as clear() does change the size back to 0
                                goto needbreak;
                            }
                        }
                    }
                    needbreak:;
                }else if(line.find("<geneRef") != std::string::npos){
                    size_t pos = line.find("id=") + 4;
                    size_t pos2 = line.find('"', pos);
                    std::string temp = line.substr(pos, pos2 - pos);
                    int indx = std::stoi(temp);
                    if(gene_id[indx - 1] == -1){ // if gene_id is -1, then it is missing data and we need to skip it
                        continue;
                    }
                    int genomeOrder{-1}; // safety measure to ensure genomeOrder exists
                    for(int i = 0; i < genomeFirstNumb.size(); ++i){
                        if(gene_id[indx - 1] >= genomeFirstNumb[i] && gene_id[indx - 1] <= genomeLastNumb[i]){
                            genomeOrder = i; // which genome the gene belongs to
                            break;
                        }
                    }
                    if(genomeOrder != -1){
                        //std::cout << "pushing back" << '\t' << genomeOrder << '\t' << geneId[indx - 1] << '\n';
                        temp_vector[genomeOrder].push_back(geneId[indx - 1]); // push on gene name at correct genomeOrder
                    }else{
                        //std::cout << "genomeOrder is -1" << '\n';
                    }
                }else if(line.find("</groups>") != std::string::npos && temp_vector.size() > 0){
                    break;
                }
            }
        }else if(line.find("</database>") != std::string::npos){
            inDB = 0;
        }
    }
    orthoxml.close();
    
    for(int i = 0; i < orthoGroup.size(); ++i){
        for(int j = 0; j < orthoGroup[i].size(); ++j){
            for(int k = 0; k < orthoGroup[i][j].size(); ++k){
                //std::cout << i << '\t' << j << '\t' << k << '\t' << orthoGroup[i][j][k] << '\n';
            }
        }
    }
    
    std::ofstream output;
    output.open(outputFileName); // Printing format to match OrthoFinder output
    output << "HOG" << '\t' << "OG" << '\t' << "Gene Tree" << '\t' << "Parent Clade" << '\t';
    for(int i = 0; i < fullInputNames.size(); ++i){
        output << fullInputNames[i] << '\t';
    }
    output << '\n';
    for(int i = 0; i < orthoGroup.size(); ++i){
        int tabCount{0}; // keeps track of how many tabs are printed a HOG line. Missing data is still printed as nothing with a tab on each side. Can add more than one "blank" tab based on the number of missing data points
        output << "N0.HOG" << std::setfill('0') << std::setw(7) << i << '\t' << "OG0000000" << '\t' << "n0" << '\t'; // Orthogroup and Node are set to 0 for all samples as we don't have that information and don't need it for OrthoRefine
        for(int j = 0; j < orthoGroup[i].size(); ++j){
            for(int k = 0; k < orthoGroup[i][j].size(); ++k){
                if(j == orthoGroup[i].size() - 1 && k == orthoGroup[i][j].size() - 1){ // if final gene from final genome
                    while(tabCount != j){
                        output << '\t';
                        ++tabCount;
                    }
                    output << orthoGroup[i][j][k]; // need this as we don't want a tab after the final gene but a newline
                }else if(orthoGroup[i][j].size() > 1 && k != orthoGroup[i][j].size() - 1){ // if there are multiple genes in a genome assigned to the same HOG, insert "," between them
                    output << orthoGroup[i][j][k] << ", ";
                }else if(tabCount != j){
                    while(tabCount != j){
                        output << '\t';
                        ++tabCount;
                    }
                    output << orthoGroup[i][j][k] << '\t';;
                    ++tabCount;
                }else{
                    output << orthoGroup[i][j][k] << '\t';
                    ++tabCount;
                }
            }
        }
        while(tabCount != fullInputNames.size() - 1){ // add additional tabs if there are missing data points after all genes that are in vector have been printed
            output << '\t';
            ++tabCount;
        }
        output << '\n';
    }
    output.close();
}