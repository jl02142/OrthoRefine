# OrthoRefine: automated enhancement of prior ortholog identification via synteny. 

[Quickstart]()\
[Install](https://github.com/jl02142/OrthoRefine/tree/main#install)\
[Input]()\
[Running]()\
[Window size & synteny ratio]()

## Quickstart
OrthoRefine may be installed on Linux systems with a C++ compiler

`
g++ -O3 orthorefine.cpp -o orthorefine.exe
`

Example required user created input file, "input.txt", of GCF accession per genome per line

>GCF_000005845.2\
>GCF_013892435.1\
>GCF_016904755.1\
>GCF_902709585.1

Single command run

`
./master_OrthoRefine.sh --input input.txt --OF_file N0.tsv --window_size window_size_number --synteny_ratio synteny_ratio_number --OrthoRefine orthorefine.exe --OrthoFinder /path/to/orthofinder.exe
`

## Install
### OrthoFinder
OrthoFinder may be found at its [Github page](https://github.com/davidemms/OrthoFinder)

### OrthoRefine
OrthoRefine's install files may be found at its [Github release page](https://github.com/jl02142/OrthoRefine/releases)

OrthoRefine may be installed on Linux systems with a C++ compiler

`
g++ -O3 orthorefine.cpp -o orthorefine.exe
`

OrthoRefine may be installed on Windows systems; an easy way to access Ubuntu on Windows is the [Windows Subsystem for Linux](https://learn.microsoft.com/en-us/windows/wsl/install). Depending on the compiler and version, additional compiler options may be required.

`
g++ -std=c++17 -pthread -O3 orthorefine.cpp -o orthorefine.exe
`

## Required input
As input, OrthoRefine requires OrthoFinder's output ("N0.tsv"), NCBI RefSeq feature table file per genome, and a user created text file where each line contains the GCF accession per genome; we reccomend having all inputs and the executable in a single directory to reduce errors in path. If a feature table file is not available for a genome, the data may be submitted to NCBI for annotation or the user may generate the annotation using the same pipeline as NCBI, [pgap](https://github.com/ncbi/pgap).

Example user created input file, "input.txt", of GCF accession per genome per line

>GCF_000005845.2\
>GCF_013892435.1\
>GCF_016904755.1\
>GCF_902709585.1

## Running OrthoRefine

OrthoRefine may be run indepdently of OrthoFinder or a support script, [master_OrthoRefine.sh](link), may be used to download the data files and run OrthoFinder and OrthoRefine with a single command. 

As single command

`
./master_OrthoRefine.sh --input input.txt --OF_file N0.tsv --window_size window_size_number --synteny_ratio synteny_ratio_number --OrthoRefine orthorefine.exe --OrthoFinder /path/to/orthofinder.exe
`

Indepedently

`
./orthorefine.exe --input input.txt --OF_file N0.tsv --window_size window_size_number --synteny_ratio synteny_ratio_number 
`

The support script, [download_for_blast.sh](link), may be run independently to download the fasta and feature tables files. 

### Runtime parameters (window size & synteny ratio)
To determine default parameters, we evaluated different combinations of window size and synteny ratio using average max number of orthologous genes (AMNOG) represented in a single SOG (syntenous ortholog group). We reccomend using a smaller window size (default 8) and higher synteny ratio (0.5), espically for closely related genomes; a larger window or lower synteny ratio may be better suited as the evolutionary distance of the genomes increases. Users may view the AMNOG for their dataset on predetermined combinations of window size and synteny ratio by not providing the window size and synteny ratio at runtime. 
