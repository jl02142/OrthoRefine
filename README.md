# OrthoRefine: automated enhancement of prior ortholog identification via synteny. 

[Quickstart](https://github.com/jl02142/OrthoRefine#quickstart)\
[OrthoRefine method summary](https://github.com/jl02142/OrthoRefine#orthorefines-method-summary)\
[Install](https://github.com/jl02142/OrthoRefine/tree/main#install)\
[Input](https://github.com/jl02142/OrthoRefine#required-input)\
[Running](https://github.com/jl02142/OrthoRefine#running-orthorefine)\
[Window size & synteny ratio](https://github.com/jl02142/OrthoRefine#runtime-parameters-window-size--synteny-ratio)

## Quickstart
OrthoRefine may be installed on Linux systems with a C++ compiler

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

## OrthoRefine's method summary
OrthoRefine automates using synteny (conserved gene order) information to refine prior ortholog identification. The analysis begins by constructing a window of user specified size centered at each gene of the HOG. Excluding this gene from the HOG (located at the center of the window), OrthoRefine evaluates the synteny by counting matching pairs of genes inside the window; matching pairs consist of genes assigned to the same HOG in the initial OrthoFinder output (Figure 1). We note that genes only need to be within the window and are not required to be in the same order, and genes that do not have a homolog in the other genome are not included in the window. The synteny ratio is calculated by taking the number of matching pairs and dividing it by the window size. If the ratio is greater than a cutoff (default 0.5), the genes at the center of the window are considered syntenic. 

<figure>
    <img src="https://github.com/jl02142/OrthoRefine/assets/23033795/9329a402-7014-4e37-909c-7531b9d45b00" width="1000" height="400">
    <figcaption>Figure 1. The window around <em>E. fergusonii’s</em> HVX45_RS11505 and its OrthoFinder matches b3271 & b0652 of <em>E. coli</em>; all three were assigned to HOG 19 by OrthoFinder and are represented by yellow circles with the red box drawn around them. Other orthologous genes (those assigned to matching HOGs by OrthoFinder) within the 10-gene window are designated by the same colored circles. In contrast, white circles indicate that the gene has no ortholog within the same window in the other genomes. The first number below each circle is the HOG assigned by OrthoFinder, while the second entry is the locus tag. As nine out of ten genes surrounding RS11505 had a match in the window centered at b3271, we concluded that there is a syntenic relationship between <em>E. coli</em> b3271 and <em>E. fergusonii</em> RS11505, and they are orthologs while b0652 is presumed to be a paralog of RS11505; none of the genes surrounding b0562 had a match within the window around RS11505 </figcaption>
</figure>

\
For additional information and examples, see OrthoRefine's paper. 

[Ludwig, J and Mrázek, J. OrthoRefine: automated enhancement of prior ortholog identification via synteny. 2023]()


## Install
### OrthoFinder
OrthoFinder may be found at its [Github page](https://github.com/davidemms/OrthoFinder)

If you use OrthoFinder, please cite their paper(s).

[Emms, D.M. and Kelly, S. (2019) OrthoFinder: phylogenetic orthology inference for comparative genomics. Genome Biology 20:238](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1832-y)

[Emms, D.M. and Kelly, S. (2015) OrthoFinder: solving fundamental biases in whole genome comparisons dramatically improves orthogroup inference accuracy. Genome Biology 16:157](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0721-2)

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
To determine default parameters, we evaluated different combinations of window size and synteny ratio using average max number of orthologous genes (AMNOG) represented in a single SOG (syntenous ortholog group). We reccomend using a smaller window size (default 8) and higher synteny ratio (default 0.5), espically for closely related genomes; a larger window or lower synteny ratio may be better suited as the evolutionary distance of the genomes increases. Users may view the AMNOG for their dataset on predetermined combinations of window size and synteny ratio by not providing the window size and synteny ratio at runtime. 
