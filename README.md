# OrthoRefine: automated enhancement of prior ortholog identification via synteny. 

[Quickstart](https://github.com/jl02142/OrthoRefine#quickstart)\
[OrthoRefine method summary](https://github.com/jl02142/OrthoRefine#orthorefines-method-summary)\
[Install](https://github.com/jl02142/OrthoRefine/tree/main#install)\
[Input](https://github.com/jl02142/OrthoRefine#required-input)\
[Running](https://github.com/jl02142/OrthoRefine#running-orthorefine)\
[Window size & synteny ratio](https://github.com/jl02142/OrthoRefine#runtime-parameters-window-size--synteny-ratio)\
[Runtime options](https://github.com/jl02142/OrthoRefine/tree/main#runtime-options)\
[Interpreting the ouput](https://github.com/jl02142/OrthoRefine/tree/main#interpreting-the-ouput)\
[Support Scripts](https://github.com/jl02142/OrthoRefine#support-scripts)

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

Single command to run OrthoFinder and OrthoRefine\
`
./master_OrthoRefine.sh --input input.txt --OF_file N0.tsv --window_size window_size_number --synteny_ratio synteny_ratio_number --OrthoRefine orthorefine.exe --OrthoFinder /path/to/orthofinder.exe -f /path/to/fasta
`

Command to run OrthoRefine only\
`
./orthorefine.exe --input input.txt --OF_file N0.tsv --window_size window_size_number --synteny_ratio synteny_ratio_number
`

The default outfile is called "outfile_ws_sr_pa_ra" where ws is the value of window size, sr is the synteny ratio, pa is print_all (default 0), and ra is run_all_orthofinder (default 0); these four suffixes are always appended. OrthoRefine's outfile has been formatted to match OrthoFinder's, except the second column will print the SOG number and the third will print the gene name from the feature table of the first genome in the SOG (genomes are ordered by their apperance in the user created input file). 

## OrthoRefine's method summary
OrthoRefine automates using synteny (conserved gene order) information to refine prior homolog (orthologous group) identification. The analysis begins by constructing a window of user specified size centered at each gene of the HOG. Excluding this gene from the HOG (located at the center of the window), OrthoRefine evaluates the synteny by counting matching pairs of genes inside the window; matching pairs consist of genes assigned to the same group by a prior program (e.g. HOG group in the OrthoFinder output). We note that genes only need to be within the window and are not required to be in the same order, and genes that do not have a homolog in the other genome are not included in the window (the window would be extended by one per missing homolog for that pairwise comparison). The synteny ratio is calculated by taking the number of matching pairs and dividing it by the window size. If the ratio is greater than a cutoff (default 0.5), the genes at the center of the window are considered syntenic. 

<figure>
    <img src="https://github.com/jl02142/OrthoRefine/assets/23033795/8f711260-18b7-4681-a5e4-52021f741206" width="1000" height="400">
    <figcaption>The window around three genes assigned to HOG19 by OrthoFinder demonstrates how OrthoRefine determines which of the E. coli genes is an ortholog of E. fergusonii’s HVX45_RS11505. The HOG19 genes are shown with yellow fill, other genes assigned to the same HOG are shown in matching colors, and genes that have orthologs in other genomes outside the displayed window are shown in white. The first number below each circle denotes the HOG assigned by OrthoFinder, while the second entry shows the locus tag.</figcaption>
</figure>

\
For additional information and examples, see OrthoRefine's paper. 

[Ludwig, J and Mrázek, J. OrthoRefine: automated enhancement of prior ortholog identification via synteny. 2023]()


## Install
### OrthoFinder
OrthoRefine requires, as input, a file of prior homolog identification; OrthoRefine currently supports OrthoFinder's default output. OrthoFinder may be found at its [Github page](https://github.com/davidemms/OrthoFinder)

### OrthoRefine
OrthoRefine's install files may be found at its [Github release page](https://github.com/jl02142/OrthoRefine/releases)

OrthoRefine can be installed in GNU or UNIX shell on Linux or MacOS

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

### Eukaryote data
While we found OrthoRefine to function on eukaryote data (<i>Saccharomyces</i>), it will not function with all eukaryote datasets. OrthoRefine currently requires the "locus_tag" column to contain data in the RefSeq feature table file. Some eukaryote data at RefSeq is missing the "locus_tag" information. Additonally, OrthoRefine currently does not handle the repeated gene identifier from isoforms in the annotation. A planned future update would resolve these issues. 

A modified input file may be submitted to tell OrthoRefine if a genome is linear or circular (second column c or l) and if it is archaea, bacteria, or eukaryote (third column a, b, or e). Circular genomes have their ends compared for syntenty (the window can "overflow" from one end to the other) while linear genes do not; by default, OrthoRefine analyzes all genomes as circular. Genomes denoted as archaea or bacteria will have operons detected by the gene gap method [(Yan and Moult. 2006.)](https://pubmed.ncbi.nlm.nih.gov/16755590/); by default, OrthoRefine does not consider operons. An operon may only count once per window for a match regardless of how many genes in the operon would have matched, the window is extended to account for this. 

>GCF_000005845.2 c b\
>GCF_013892435.1 c b\
>GCF_016904755.1 c b\
>GCF_902709585.1 l e

## Running OrthoRefine

OrthoRefine may be run indepdently of OrthoFinder or a support script, [master_OrthoRefine.sh](https://github.com/jl02142/OrthoRefine/blob/main/master_OrthoRefine.sh), may be used to download the data files and run OrthoFinder and OrthoRefine with a single command. 

As single command

`
./master_OrthoRefine.sh --input input.txt --OF_file N0.tsv --window_size window_size_number --synteny_ratio synteny_ratio_number --OrthoRefine orthorefine.exe --OrthoFinder /path/to/orthofinder.exe -f /path/to/fasta
`

Indepedently

`
./orthorefine.exe --input input.txt --OF_file N0.tsv --window_size window_size_number --synteny_ratio synteny_ratio_number 
`

### Runtime parameters (window size & synteny ratio)
To determine default parameters, we evaluated different combinations of window size and synteny ratio using average max number of orthologous genes (AMNOG) represented in a single SOG (syntenous ortholog group). We reccomend using a smaller window size (default 8) and higher synteny ratio (default 0.5), espically for closely related genomes; a larger window or lower synteny ratio may be better suited as the evolutionary distance of the genomes increases (e.g., window size 30 and synteny ratio 0.2). Users may view the AMNOG for their dataset on predetermined combinations of window size and synteny ratio by setting the runtime option --run_combo to 1 (Calculating the AMNOG is a parallelized process that is currently memory intensive. OrthoRefine will return the killed error if the memory of the system is exceeded). 

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

--diag will print extra information to diagnosis potential problems. Accepted values are 0 (none, default), 1 (short), or 2 (long) or 3 (long long) or 4 (everything). Setting diag to 2 will allow users to see which genes were in the window and which of those were matches; we recommend that users combine --diag 2 with --run_single_HOG to reduce the amount of text printed to the std out. 

### Interpreting the ouput

OrthoRefine's output has been formatted to closely match OrthoFidner's. A change has been made to the second collumn where the orthogroup has been replaced with the SOG, and to the third collumn where the node has been replaced with the gene name from the feature table file of the first genome listed in the input file. The final line of the output contains the number of HOGs refined and the total number of refinements, which can be larger as a single HOG can be refined into mulitple SOGs. 

An example of OrthoRefine's output is below. HOG 5 is split into two SOGs, 5.0 which contains b1552 and HVX45_RS19450 and 5.1 which contains b0990, JRC41_RS13205, and GV529_RS02480.
```
HOG     SOG     Gene_name       GCF_000005845.2_ASM584v2_feature_table.txt      GCF_013892435.1_ASM1389243v1_feature_table.txt  GCF_016904755.1_ASM1690475v2_feature_table.txt  GCF_902709585.1_H1-003-0086-C-F.v2_feature_table.txt
N0.HOG0000005   5.0     cold shock-like protein CspI    b1552   HVX45_RS19450
N0.HOG0000005   5.1     cold shock protein CspG b0990           JRC41_RS13205   GV529_RS02480
N0.HOG0000019   19.0    glutamate/aspartate ABC transporter ATP binding subunit b0652   HVX45_RS07420   JRC41_RS15115   GV529_RS05870
N0.HOG0000019   19.1    putative ABC transporter ATP-binding subunit YhdZ       b3271   HVX45_RS11505           GV529_RS14465
...
Number of HOGs refined: 408     for a total refinement of       468
```

## Support Scripts

There are 3 Bash support scripts to support OrthoRefine: [master_OrthoRefine.sh](https://github.com/jl02142/OrthoRefine/blob/main/master_OrthoRefine.sh) may be used to download the data files and run OrthoFinder and OrthoRefine with a single command, [download_ft_fafiles.sh](https://github.com/jl02142/OrthoRefine/blob/main/download_ft_fafiles.sh) may be run independently to download the fasta and feature tables files, and [summary_stats.sh](https://github.com/jl02142/OrthoRefine/blob/main/summary_stats.sh) may be used to generate summary stats (number of 1-to-1 HOGs, 0-or-1 HOGS, paralogous HOGs. number of HOGs confirmed, unconfirmed, or modified by synteny) for a particular window size and synteny ratio. 

`
./master_OrthoRefine.sh --input input.txt --OF_file N0.tsv --window_size window_size_number --synteny_ratio synteny_ratio_number --OrthoRefine orthorefine.exe --OrthoFinder /path/to/orthofinder.exe -f /path/to/fasta
`

`
./download_ft_fafiles.sh input.txt
`

`
./summary_stats.sh --input input.txt --OF_file N0.tsv --window_size 8 --synteny_ratio 0.5 --exe ./OrthoRefine.exe
`

The Linux commands chmod or dos2unix may be required to use the support scripts. chmod changes the script permissions so it may be run. dos2unix will make sure the file is in Linux format and not Windows (dos) format. 

`
chmod u+x download_ft_fafiles.sh
`

`
dos2unix download_ft_fafiles.sh
`
