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

The default outfile is called "outfile_ws_sr_pa_ra" where ws is the value of window size (default 8), sr is the synteny ratio (default 0.5), pa is print_all (default 0), and ra is run_all_orthofinder (default 0); these four suffixes are always appended and the default output file name is "outfile_8_0.5_0_0". OrthoRefine's outfile has been formatted to match OrthoFinder's, except the second column will print the SOG number and the third will print the gene name from the feature table of the first genome in the SOG (genomes are ordered by their apperance in the user created input file). 

## OrthoRefine's method summary
OrthoRefine is a tool designed to automate the refinement of hierarchical orthogroups (HOGs) identification using synteny, conservation of gene order. The analysis proceeds through the following steps:

* Construction of Gene Windows: OrthoRefine begins by constructing a window of user-specified size centered at each gene of the HOG. This window excludes the gene located at the center and evaluates synteny by counting matching pairs of genes within the window.

* Evaluation of Synteny: Matching pairs consist of genes assigned to the same group by a prior program (e.g., HOG group in the OrthoFinder output). It's important to note that genes only need to be within the window and are not required to be in the same order. Additionally, genes without a homolog in the other genome are not included in the window, and the window is extended by one per missing homolog for that pairwise comparison.

* Calculation of Synteny Ratio: The synteny ratio is calculated by dividing the number of matching pairs by the window size. If the ratio exceeds a user-defined cutoff (default 0.5), the genes at the center of the window are considered syntenic.

OrthoRefine provides a refined assessment of homologs based on synteny information, which can aid in improving the accuracy of orthologous group assignments. Users can adjust parameters such as the window size and synteny ratio cutoff to customize the analysis according to their specific requirements.

<figure>
    <img src="https://github.com/jl02142/OrthoRefine/assets/23033795/8f711260-18b7-4681-a5e4-52021f741206" width="1000" height="400">
    <figcaption>The window around three genes assigned to HOG19 by OrthoFinder demonstrates how OrthoRefine determines which of the <i>E. coli</i> genes is an ortholog of <i>E. fergusonii’s</i> HVX45_RS11505. The HOG19 genes are shown with yellow fill, other genes assigned to the same HOG are shown in matching colors, and genes that have orthologs in other genomes outside the displayed window are shown in white. The first number below each circle denotes the HOG assigned by OrthoFinder, while the second entry shows the locus tag.</figcaption>
</figure>

\
OrthoRefine's manuscript is available at BMC Bioinformatics. 

[Ludwig, J and Mrázek, J. OrthoRefine: automated enhancement of prior ortholog identification via synteny.](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-024-05786-7)


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
As input, OrthoRefine requires OrthoFinder's output ("N0.tsv") which contains orthogroup information, NCBI RefSeq feature table file per genome which provides the genome annotation information, and a user created text file where each line contains the GCF accession per genome. We reccomend having all inputs and the executable in a single directory to reduce errors in path. If a feature table file is not available for a genome, the data may be submitted to NCBI for annotation or the user may generate the annotation using the same pipeline as NCBI, [pgap](https://github.com/ncbi/pgap). OrthoFinder's output (N0.tsv) is in Window's format (\r\n) - to do the necesarry conversion to Linux (\n), use the Linux command dos2unix; the master script already handles this conversion.

`
dos2unix N0.tsv
`

Example user created input file, "input.txt". Each line must contain one GCF accession.

>GCF_000005845.2\
>GCF_013892435.1\
>GCF_016904755.1\
>GCF_902709585.1

### Eukaryote data
While OrthoRefine has been successfully tested on eukaryotic datasets such as <i>Saccharomyces</i>, it may encounter limitations with certain eukaryote datasets. Currently, OrthoRefine requires the 'locus_tag' column to contain data in the RefSeq feature table file. However, some eukaryotic data available at RefSeq may lack 'locus_tag' information, which can prevent OrthoRefine from functioning properly.

Additionally, OrthoRefine does not handle repeated gene identifiers from isoforms in the annotation files. To address this issue, users should utilize the provided Python script [ft_fa_isoform_remove.py](https://github.com/jl02142/OrthoRefine/blob/main/ft_fa_isoform_remove.py) to remove isoforms from both the fasta and feature table files before initiating the analysis.

#### Circular vs linear genomes and operon detection

A modified input file format may be submitted to provide additional information to OrthoRefine, allowing users to specify if a genome is linear or circular and if it belongs to the archaea, bacteria, or eukaryote domain. This information is provided in the second and third columns of the input file:

* The second column denotes the genome's topology, with 'c' indicating a circular genome and 'l' indicating a linear genome. Circular genomes have their ends compared for synteny, meaning that the window can 'overflow' from one end to the other. By default, OrthoRefine analyzes all genomes as circular unless otherwise specified.

* The third column specifies the domain of the genome, with 'a' representing archaea, 'b' representing bacteria, and 'e' representing eukaryote.

Genomes denoted as archaea or bacteria will have operons detected using the gene gap method described by [Yan and Moult (2006)](https://pubmed.ncbi.nlm.nih.gov/16755590/). By default, OrthoRefine does not consider operons. Operons are counted once per window for a match, regardless of the number of genes in the operon that would have matched. The window is extended to account for this.

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
To establish default parameters, we conducted an evaluation of various combinations of window size and synteny ratio, assessing the Average Max Number of Orthologous Genes (AMNOG) represented in a single Syntenous Ortholog Group (SOG). Based on our analysis, we recommend defaulting to a smaller window size (set at 8) and a higher synteny ratio (set at 0.5), particularly for closely related genomes. However, as the evolutionary distance between genomes increases, users may find that a larger window size or a lower synteny ratio is more appropriate (e.g., a window size of 30 and a synteny ratio of 0.2).

Users have the option to examine the AMNOG for their dataset across predetermined combinations of window size and synteny ratio by enabling the runtime option --run_combo to 1. It's important to note that calculating the AMNOG is a parallelized process that currently requires significant memory resources. If the system's memory capacity is exceeded during this calculation, OrthoRefine will terminate with a 'killed' error.

In addition to adjusting the window size and synteny ratio parameters, users must specify the input file (--input) containing user-created datafile of GCF accession per genome and the output file from OrthoFinder (--OF_file) to initiate the analysis.

### Runtime options
OrthoRefine has several runtime options, some which standard end-users may find useful and others intended for advanced end-users. 

#### Standard options
```
--print_all
```
Controls if OrthoRefine should print only those groups with changes supported by synteny (0, default) or print all groups even if no synteny support (1), or print only groups supported by synteny - even if no change occured in the group (2). \
```
--run_all_orthofinder
```
Controls if OrthoRefine should only evaluate synteny on groups with paralogs (0, default) or on all groups (1).\
```
--outfile
```
Sets the prefix of the outfile. The default outfile is called "outfile_ws_sr_pa_ra" where ws is the value of window size, sr is the synteny ratio, pa is print_all (default 0), and ra is run_all_orthofinder (default 0); these four suffixes are always appended. \
```
--run_single_HOG
```
 Controls if OrthoRefine should only evaluate a single group for synteny, specified by the value. \
```
--prod_acc
```
Controls if OrthoRefine should print the product accession (0) or locus tag (1, default). We prefer the locus tag over product accession as the locus tag is non-redundant while the product accession may be redundant.\
```
--path
```
File path to location if the input files are located in a different location than the executable.

#### Advanced options

--diag will print extra information to diagnosis potential problems. Accepted values are 0 (none, default), 1 (short), or 2 (long) or 3 (long long) or 4 (everything). Setting diag to 2 will allow users to see which genes were in the window and which of those were matches; we recommend that users combine --diag 2 with --run_single_HOG to reduce the amount of text printed to the std out. 

### Interpreting the ouput

OrthoRefine's output closely resembles that of OrthoFinder, with some modifications. Notably, the second column now displays the Syntenous Ortholog Group (SOG) instead of the orthogroup, while the third column contains the gene name from the feature table file of the first genome listed in the input file, replacing the node identifier.

The final line of the output provides summary statistics, including the number of Homologous Orthologous Groups (HOGs) refined and the total number of refinements. It's important to note that the total number of refinements may exceed the number of refined HOGs, as a single HOG can be refined into multiple SOGs.

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
