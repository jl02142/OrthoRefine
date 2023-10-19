# OrthoRefine: automated enhancement of prior ortholog identification via synteny. 

[Quickstart](https://github.com/jl02142/OrthoRefine#quickstart)\
[OrthoRefine method summary](https://github.com/jl02142/OrthoRefine#orthorefines-method-summary)\
[Install](https://github.com/jl02142/OrthoRefine/tree/main#install)\
[Input](https://github.com/jl02142/OrthoRefine#required-input)\
[Running](https://github.com/jl02142/OrthoRefine#running-orthorefine)\
[Window size & synteny ratio](https://github.com/jl02142/OrthoRefine#runtime-parameters-window-size--synteny-ratio)\
[Runtime options](https://github.com/jl02142/OrthoRefine/tree/main#runtime-options)

## Quickstart
OrthoRefine may be installed with a C++ compiler

`
g++ -O3 orthorefine.cpp -o orthorefine.exe
`

Example required user created input file, "input.txt". Each line must contain one GCF accession.

>GCF_000005845.2\
>GCF_013892435.1\
>GCF_016904755.1\
>GCF_902709585.1

Single command run

`
./master_OrthoRefine.sh --input input.txt --OF_file N0.tsv --window_size window_size_number --synteny_ratio synteny_ratio_number --OrthoRefine orthorefine.exe --OrthoFinder /path/to/orthofinder.exe
`

The default outfile is called "outfile_ws_sr_pa_ra" where ws is the value of window size, sr is the synteny ratio, pa is print_all (default 0), and ra is run_all_orthofinder (default 0); these four suffixes are always appended. 

## OrthoRefine's method summary
OrthoRefine automates using synteny (conserved gene order) information to refine prior homolog (orthologous group) identification. The analysis begins by constructing a window of user specified size centered at each gene of the HOG. Excluding this gene from the HOG (located at the center of the window), OrthoRefine evaluates the synteny by counting matching pairs of genes inside the window; matching pairs consist of genes assigned to the same group by a prior program (e.g. HOG group in the OrthoFinder output) (Figure 1). We note that genes only need to be within the window and are not required to be in the same order, and genes that do not have a homolog in the other genome are not included in the window. The synteny ratio is calculated by taking the number of matching pairs and dividing it by the window size. If the ratio is greater than a cutoff (default 0.5), the genes at the center of the window are considered syntenic. 

<figure>
    <img src="https://github.com/jl02142/OrthoRefine/assets/23033795/9329a402-7014-4e37-909c-7531b9d45b00" width="1000" height="400">
    <figcaption>Figure 1. The window around <em>E. fergusonii’s</em> HVX45_RS11505 and its OrthoFinder matches b3271 & b0652 of <em>E. coli</em>; all three were assigned to HOG 19 by OrthoFinder and are represented by yellow circles with the red box drawn around them. Other orthologous genes (those assigned to matching HOGs by OrthoFinder) within the 10-gene window are designated by the same colored circles. In contrast, white circles indicate that the gene has no ortholog within the same window in the other genomes. The first number below each circle is the HOG assigned by OrthoFinder, while the second entry is the locus tag. As nine out of ten genes surrounding RS11505 had a match in the window centered at b3271, we concluded that there is a syntenic relationship between <em>E. coli</em> b3271 and <em>E. fergusonii</em> RS11505, and they are orthologs while b0652 is presumed to be a paralog of RS11505; none of the genes surrounding b0562 had a match within the window around RS11505 </figcaption>
</figure>

\
For additional information and examples, see OrthoRefine's paper. 

[Ludwig, J and Mrázek, J. OrthoRefine: automated enhancement of prior ortholog identification via synteny. 2023]()


## Install
### OrthoFinder
OrthoRefine requires, as input, a file of prior homolog identification; OrthoRefine currently supports OrthoFinder's default output. OrthoFinder may be found at its [Github page](https://github.com/davidemms/OrthoFinder)

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
As input, OrthoRefine requires OrthoFinder's output ("N0.tsv"), NCBI RefSeq feature table file per genome, and a user created text file where each line contains the GCF accession per genome; we reccomend having all inputs and the executable in a single directory to reduce errors in path. If a feature table file is not available for a genome, the data may be submitted to NCBI for annotation or the user may generate the annotation using the same pipeline as NCBI, [pgap](https://github.com/ncbi/pgap). OrthoFinder's output (N0.tsv) is in Window's format (\r\n) - to do the necesarry conversion to Linux (\n), use the Linux command dos2unix; the master script already handles this.

`
dos2unix N0.tsv
`

Example user created input file, "input.txt". Each line must contain one GCF accession.

>GCF_000005845.2\
>GCF_013892435.1\
>GCF_016904755.1\
>GCF_902709585.1

## Running OrthoRefine

OrthoRefine may be run indepdently of OrthoFinder or a support script, [master_OrthoRefine.sh](https://github.com/jl02142/OrthoRefine/blob/main/master_OrthoRefine.sh), may be used to download the data files and run OrthoFinder and OrthoRefine with a single command. 

As single command

`
./master_OrthoRefine.sh --input input.txt --OF_file N0.tsv --window_size window_size_number --synteny_ratio synteny_ratio_number --OrthoRefine orthorefine.exe --OrthoFinder /path/to/orthofinder.exe
`

Indepedently

`
./orthorefine.exe --input input.txt --OF_file N0.tsv --window_size window_size_number --synteny_ratio synteny_ratio_number 
`

The support script, [download_ft_fafiles.sh](https://github.com/jl02142/OrthoRefine/blob/main/download_ft_fafiles.sh), may be run independently to download the fasta and feature tables files. 

### Runtime parameters (window size & synteny ratio)
To determine default parameters, we evaluated different combinations of window size and synteny ratio using average max number of orthologous genes (AMNOG) represented in a single SOG (syntenous ortholog group). We reccomend using a smaller window size (default 8) and higher synteny ratio (default 0.5), espically for closely related genomes; a larger window or lower synteny ratio may be better suited as the evolutionary distance of the genomes increases. Users may view the AMNOG for their dataset on predetermined combinations of window size and synteny ratio by not providing the window size and synteny ratio at runtime. 

Additional runtime parameters are --input, the user created input file, and --OF_file, the output from OrthoFinder. 

### Runtime options
OrthoRefine has several runtime options, some which standard end-users may find useful and others intended for advanced end-users. 

#### Standard options
--print_all Controls if OrthoRefine should print only those groups with changes supported by synteny (0, default) or print all groups even if no synteny support (1), or print only groups supported by synteny - even if no change occured in the group (2). 
--run_all_orthofinder Controls if OrthoRefine should only evaluate synteny on groups with paralogs (0, default) or on all groups (1).
--outfile sets the prefix of the outfile. The default outfile is called "outfile_ws_sr_pa_ra" where ws is the value of window size, sr is the synteny ratio, pa is print_all (default 0), and ra is run_all_orthofinder (default 0); these four suffixes are always appended. 
--run_single_HOG Controls if OrthoRefine should only evaluate a single group for synteny, specified by the value. 
--prod_acc Controls if OrthoRefine should print the product accession (0, default) or locus tag (1).
--path File path to location if the input files are located in a different location than the executable.

#### Advanced options
--diag will print extra information to diagnosis potential problems. Accepted values are 0 (none, default), 1 (short), or 2 (long) or 3 (long long) or 4 (everything). Not intended for everyday use. 
