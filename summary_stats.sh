# from https://www.codebyamir.com/blog/parse-command-line-arguments-using-getopt
window_size=8
synteny_ratio=0.5
SHORT=:
LONG=input:,OF_file:,window_size:,synteny_ratio:,exe:,sankey_out:,core:,complex:,
OPTS=$(getopt --options $SHORT --long $LONG --name "$0" -- "$@")

if [ $? != 0 ]; then echo "ERROR: bad option provided." >&2 ; exit 1; fi

eval set -- "$OPTS"

complex=0

while true; do
    case "$1" in
        - | --input )
            input="$2" # input file with REFSEQ GCF for OrthoRefine
            shift 2
            ;;
        - | --OF_file ) # output file from OrthoFinder
            OF_file="$2"
            shift 2
            ;;
        - | --window_size )
            window_size="$2"
            shift 2
            ;;
        - | --synteny_ratio )
            synteny_ratio="$2"
            shift 2
            ;;
        - | --exe ) # name of OrthoRefine compiled program
            exe="$2"
            shift 2
            ;;
        - | --sankey_out ) # if user wants a different outfile name of the figure commands
            sankey_out="$2"
            shift 2
            ;;
        - | --core ) # number of cores
            core="$2"
            shift 2
            ;;
        - | --complex ) # if the more complex figure output should be printed
            complex="$2"
            shift 2
            ;;
        -- )
            shift
            break
            ;;
        * )
            echo "Error"
            exit
            ;;
    esac
done

# require user to provide OrthoRefine exe name
if [ -z "$exe" ]; then 
    echo "Error: exe not provided"
    return 1 2>/dev/null
    { exit 1; }
fi

if [ -z "$OF_file" ];then
    echo "Error: OF_file (N0.tsv) not provided"
    return 1 2>/dev/null
    { exit 1; }
fi

if [ -z "$input" ];then
    echo "Error: input not provided"
    return 1 2>/dev/null
    { exit 1; }
fi

# set number of cores if user didn't provide
if [ -z "$core"]; then
    core=$(grep -c ^processor /proc/cpuinfo)
fi

#only need 6 cores max
if [ $core -gt 6 ]; then 
    core=6
fi

# if user didn't specify sankeymatic file name
if [ -z $sankey_out ]; then
    sankey_out=sankeymatic_fig.txt
fi

# don't overwrite sankeymatic files
if (( $complex > 0 )); then
    if test -f "$sankey_out"; then
        echo "$sankey_out exists. Change current file name or delete it. Exiting"
        return 1 2>/dev/null
        { exit 1; }
    fi
fi
# create array containing the commands to be run
cluster_submit=( "${cluster_submit[@]}" "./$exe --input $input --OF_file $OF_file --window_size $window_size --synteny_ratio $synteny_ratio --print_all 0 --run_all_orthofinder 0 " )
cluster_submit=( "${cluster_submit[@]}" "./$exe --input $input --OF_file $OF_file --window_size $window_size --synteny_ratio $synteny_ratio --print_all 1 --run_all_orthofinder 0 " )
cluster_submit=( "${cluster_submit[@]}" "./$exe --input $input --OF_file $OF_file --window_size $window_size --synteny_ratio $synteny_ratio --print_all 0 --run_all_orthofinder 1 " )
cluster_submit=( "${cluster_submit[@]}" "./$exe --input $input --OF_file $OF_file --window_size $window_size --synteny_ratio $synteny_ratio --print_all 1 --run_all_orthofinder 1 " )
cluster_submit=( "${cluster_submit[@]}" "./$exe --input $input --OF_file $OF_file --window_size $window_size --synteny_ratio $synteny_ratio --print_all 2 --run_all_orthofinder 0 " )
cluster_submit=( "${cluster_submit[@]}" "./$exe --input $input --OF_file $OF_file --window_size $window_size --synteny_ratio $synteny_ratio --print_all 2 --run_all_orthofinder 1 " )

# runs the above 6 commands in parallel based on number of cores
printf "%s\n" "${cluster_submit[@]}" | xargs -P $core -I {} sh -c - {}


# perfix of outfile created by OrthoRefine
outfile=outfile_"$window_size"_"$synteny_ratio"
#echo $outfile

# total HOGs from OrthoFinder
total_hogs=$(cat $OF_file | tail -n -1 | cut -f1 | sed 's/N0.HOG//'  | sed 's/0\{1,6\}//')
((total_hogs+=1))  # Because HOG count starts with index 0
#echo $total_hogs

# number of HOGs that are in a 1-to1 relationship from OrthoFinder
one_to_one_orthologs=$(cat $OF_file | awk -F '\t' '{for(i=3; i<=NF; ++i){if(i == 3){count=0}{if($i ~ /^[[:alpha:]][[:alpha:]]_[[:digit:]]+\.[[:digit:]]$/){++count}}}}{if(count ~ NF - 3){print $0}}' | wc -l)
#echo $one_to_one_orthologs

# number of HOGS with at least one paralog from OrthoFinder
hogs_with_paralogs=$(cat $OF_file | grep -c ",")
#echo $hogs_with_paralogs

# number of HOGs where each genome contributed 0 or 1 gene HOGs. No paralogs and not in 1-to-1 relationship
zero_or_one_hogs=$(($total_hogs-$one_to_one_orthologs-$hogs_with_paralogs))
#echo $zero_or_one_hogs

# some HOGs only contain 1 genome, OrthoRefine (synteny) cannot be used with only 1 genome
number_hogs_could_be_refined=$(cat $outfile\_1\_1 | tail -n -1 | cut -f2)
#echo $number_hogs_could_be_refined

#cat $outfile\_0\_0 | head -n -1 | cut -f2 | sort | uniq > zero_zero_uniq 

# run OrthoRefine printing all paralogous HOGs that can be refined
cat $outfile\_1\_0 | head -n -1 | tail +2 | cut -f1 | cut -c 7- | sed 's/0\{1,6\}//' | sort | uniq > one_zero_uniq
# run OrthoRefine printing all paralogous HOGs that were confirmed with synteny
cat $outfile\_2\_0 | head -n -1 | tail +2 | cut -f1 | cut -c 7- | sed 's/0\{1,6\}//' | sort | uniq > two_zero_uniq
# print to file_2934 those HOGs (their numbers) that are uniq of the paralogs that could be refined (not found in paralog confirmed with synteny) and those HOGs of paralogs that were also found in paralogs confirmed with synteny
# we supressed those HOGs only found in the second file
comm -2 one_zero_uniq two_zero_uniq > file_2934
# remove tab character
less file_2934 | tr -d '\t' > temptemptemp
mv temptemptemp file_2934
# from OrthoFinder output, get lines with paralogs and then the first column with the HOG number and strip the leading 0's from the HOG number
cat $OF_file | grep "," | cut -f1 | sed 's/N0.HOG//'  | sed 's/0\{1,6\}//' | sort | uniq > file_3502
# this adjustment accounts for HOGs that contain a product accession (from OrthoRefine) that is also the product accession of a different gene located in the same genome. Product accessions are not unique, unlike locus tags
# adjustment is necesary to make the math work because of the way OrthoRefine uses locus_tags instead of product accessions
adjust=$(comm -23 file_2934 file_3502 | wc -l)
#echo "adjust"$adjust

# number of HOGs that are paralogs with one genome
number_hogs_only_1_genome_no_refine_temp=$(cat $OF_file | grep "," | cut -f4- | grep -v $'[[:alnum:]]\+.*\t[[:alnum:]]\+' | wc -l)
# this is the count of paralogous HOGs that could not be refined due to only 1 genome with the adjustment
number_hogs_only_1_genome_no_refine=$((number_hogs_only_1_genome_no_refine_temp-$adjust))
#echo $number_hogs_only_1_genome_no_refine_temp
#echo $number_hogs_only_1_genome_no_refine


# number of paralogs refined with OrthoRefine
number_hogs_with_paralogs_refined=$(cat $outfile\_0\_0 | tail -n -1 | cut -f2)
#echo $number_hogs_with_paralogs_refined

#number of paralogs that could not be refined
number_hogs_with_paralogs_unrefined=$(($hogs_with_paralogs-$number_hogs_only_1_genome_no_refine-$number_hogs_with_paralogs_refined))
#echo $number_hogs_with_paralogs_unrefined

# HOGs in 1-to-1 relationship. Skipping first 3 columns as they contain no product accessions. Those in 1-to-1 will have a count equal to the number of fields - 3
cat $OF_file| awk -F '\t' '{for(i=3; i<=NF; ++i){if(i == 3){count=0}{if($i ~ /^[[:alpha:]][[:alpha:]]_[[:digit:]]+\.[[:digit:]]$/){++count}}}}{if(count ~ NF - 3){print $0}}' > for_grep_for_figure_1_to_1
# get HOG number and strip leading 0's
cat for_grep_for_figure_1_to_1 | cut -f1 | sed 's/N0.HOG//'  | sed 's/0\{1,6\}//' | sort > to_compare_1
# run OrthoRefine printing those refined by synteny and on all HOGs
cat $outfile\_0\_1 | head -n -1 | tail +2 | cut -f1 | cut -c 7- | sed 's/0\{1,6\}//' | sort | uniq > to_compare_2
# number of HOGs that were 1-1 that could be refined by OrthoRefine. Print those HOGs present in both files
number_hogs_one_to_one_refined=$(comm -12 to_compare_1 to_compare_2 | wc -l)
#echo $number_hogs_one_to_one_refined

# number of HOGs in 1-to-1 unrefined
number_hogs_one_to_one_unrefined=$(($one_to_one_orthologs-$number_hogs_one_to_one_refined))
#echo $number_hogs_one_to_one_unrefined

# number of HOGs refined by synteny
number_of_hogs_changed=$(cat $outfile\_0\_1 | tail -n -1 | cut -f2) 
#echo $number_of_hogs_changed

# number of paralogous HOGs changed with synteny (CWS)
number_of_hogs_with_paralogs_CWS=$(cat $outfile\_0\_0 | tail -n -1 | cut -f2)
#echo $number_of_hogs_with_paralogs_CWS

#number of HOGs that were 0 or 1 gene HOGs per genome refined
number_hogs_zero_or_one_refined=$(($number_of_hogs_changed-$number_of_hogs_with_paralogs_CWS-$number_hogs_one_to_one_refined))
#echo $number_hogs_zero_or_one_refined

#number of HOGs that were 0 or 1 gene HOGs per genome unrefined
number_of_hogs_zero_or_one_unrefined=$(($zero_or_one_hogs-$number_hogs_zero_or_one_refined))
#echo $number_of_hogs_zero_or_one_unrefined

# outfile_2_1 is from OrthRefine run on all HOGs and printing those HOGs that were changed with synteny
cat $outfile\_2\_1 | head -n -1 | tail +2 | cut -f1 | cut -c 7- | sed 's/0\{1,6\}//' | sort | uniq > to_compare_3
#number of HOGs that were 1-to-1 that were confirmed by synteny (Unchanged Synteny confirmed and changed with synteny)
number_of_hogs_one_to_one_orthologs_USC_and_CWS=$(comm -12 to_compare_1 to_compare_3 | wc -l)
#echo $number_of_hogs_one_to_one_orthologs_USC_and_CWS

# number of 1-to-1 HOGs unchanged and synteny unconfirmed (USU)
number_of_hogs_one_to_one_USU=$(($one_to_one_orthologs-$number_of_hogs_one_to_one_orthologs_USC_and_CWS))
#echo $number_of_hogs_one_to_one_USU

# number of 1-to-1 HOGs Unchanged confirmed with synteny (USC)
number_of_hogs_one_to_one_USC=$(($one_to_one_orthologs-$number_hogs_one_to_one_refined-$number_of_hogs_one_to_one_USU))
#echo $number_of_hogs_one_to_one_USC

# print all paralogous HOGs and print all paralgous HOGs confirmed with synteny (both changed and unchanged)
temp1=$(cat $outfile\_1\_0 | tail -n -1 | cut -f2)
temp2=$(cat $outfile\_2\_0 | tail -n -1 | cut -f2)
# number of paralogs unchanged synteny unconfirmed (USU)
number_paralogs_USU=$(($temp1-$temp2))
#echo $number_paralogs_USU

# number of paralgous unchanged but confirmed by synteny (USC)
number_paralogs_USC=$(($temp2-$number_hogs_with_paralogs_refined))
#echo $number_paralogs_USC

temp3=$(cat $outfile\_2\_1 | tail -n -1 | cut -f2)
temp4=$(($number_of_hogs_one_to_one_USC+$number_paralogs_USC+$number_hogs_one_to_one_refined+$number_hogs_with_paralogs_refined+$number_hogs_zero_or_one_refined))
number_zero_or_one_USC=$(($temp3-$temp4))
#echo $number_zero_or_one_USC

number_zero_or_one_USU=$(($number_of_hogs_zero_or_one_unrefined-$number_zero_or_one_USC))
#echo $number_zero_or_one_USU

number_hogs_USC=$(($number_of_hogs_one_to_one_USC+$number_zero_or_one_USC+$number_paralogs_USC))
#echo $number_hogs_USC

number_hogs_USU=$(($number_of_hogs_one_to_one_USU+$number_zero_or_one_USU+$number_paralogs_USU))
#echo $number_hogs_USU

number_hogs_unchanged=$(($number_hogs_USC+$number_hogs_USU))
#echo $number_hogs_unchanged


rm for_grep_for_figure_1_to_1
rm to_compare_1
rm to_compare_2
rm to_compare_3
rm file_2934
rm file_3502
rm one_zero_uniq
rm two_zero_uniq

if (( $complex > 0 )); then

    echo "OrthoFinder HOGs ["$one_to_one_orthologs"] 1-to-1 HOGs" >> $sankey_out
    echo "OrthoFinder HOGs ["$zero_or_one_hogs"] 0 or 1 gene HOGs" >> $sankey_out
    echo "OrthoFinder HOGs ["$hogs_with_paralogs"] Paralog HOGs" >> $sankey_out

    if [ $complex != 1 ]; then
        echo "Paralog HOGs ["$number_hogs_only_1_genome_no_refine"] HOGs only 1 genome #000000" >> $sankey_out
    fi

    if [ $complex == 1 ]; then 
        number_paralog_HOGs_unchanged=$(($number_hogs_with_paralogs_unrefined+$number_hogs_only_1_genome_no_refine))
        echo "Paralog HOGs ["$number_paralog_HOGs_unchanged"] HOGs unchanged #228B22" >> $sankey_out
        echo "1-to-1 HOGs ["$number_hogs_one_to_one_unrefined"] HOGs unchanged" >> $sankey_out
        echo "0 or 1 gene HOGs["$number_of_hogs_zero_or_one_unrefined"] HOGs unchanged" >> $sankey_out
        echo "Paralog HOGs ["$number_hogs_with_paralogs_refined"] HOGs changed #e58e73" >> $sankey_out
        echo "1-to-1 HOGs ["$number_hogs_one_to_one_refined"] HOGs changed #e58e73" >> $sankey_out
        echo "0 or 1 gene HOGs["$number_hogs_zero_or_one_refined"] HOGs changed #e58e73" >> $sankey_out
    elif [ $complex == 2 ]; then
        echo "Paralog HOGs ["$number_hogs_with_paralogs_refined"] HOGs changed #800080" >> $sankey_out
        echo "Paralog HOGs ["$number_paralogs_USC"] HOGs USC #800080" >> $sankey_out
        echo "Paralog HOGs ["$number_paralogs_USU"] HOGs USU #800080" >> $sankey_out
        echo "1-to-1 HOGs ["$number_hogs_one_to_one_refined"] HOGs changed" >> $sankey_out
        echo "1-to-1 HOGs ["$number_of_hogs_one_to_one_USC"] HOGs USC" >> $sankey_out
        echo "1-to-1 HOGs ["$number_of_hogs_one_to_one_USU"] HOGs USU" >> $sankey_outHOGs 0 or 1
        echo "0 or 1 gene HOGs["$number_hogs_zero_or_one_refined"] HOGs changed" >> $sankey_out
        echo "0 or 1 gene HOGs["$number_zero_or_one_USC"] HOGs USC" >> $sankey_out
        echo "0 or 1 gene HOGs["$number_zero_or_one_USU"] HOGs USU" >> $sankey_out
    elif [ $complex == 3 ]; then
        echo "Paralog HOGs ["$number_hogs_with_paralogs_refined"] Paralog HOGs CWS #800080" >> $sankey_out
        echo "Paralog HOGs ["$number_paralogs_USC"] Paralog HOGs USC #800080" >> $sankey_out
        echo "Paralog HOGs ["$number_paralogs_USU"] Paralog HOGs USU #800080" >> $sankey_out
        echo "1-to-1 HOGs ["$number_hogs_one_to_one_refined"] 1-to-1 HOGs CWS" >> $sankey_out
        echo "1-to-1 HOGs ["$number_of_hogs_one_to_one_USC"] 1-to-1 HOGs USC" >> $sankey_out
        echo "1-to-1 HOGs ["$number_of_hogs_one_to_one_USU"] 1-to-1 HOGs USU" >> $sankey_out
        echo "0 or 1 gene HOGs["$number_hogs_zero_or_one_refined"] 0 or 1 gene HOGs CWS" >> $sankey_out
        echo "0 or 1 gene HOGs["$number_zero_or_one_USC"] 0 or 1 gene HOGs USC" >> $sankey_out
        echo "0 or 1 gene HOGs["$number_zero_or_one_USU"] 0 or 1 gene HOGs USU" >> $sankey_out
        echo "Paralog HOGs CWS["$number_hogs_with_paralogs_refined"] HOGs changed #e58e73" >> $sankey_out
        echo "1-to-1 HOGs CWS ["$number_hogs_one_to_one_refined"] HOGs changed #e58e73" >> $sankey_out
        echo "0 or 1 gene HOGs CWS["$number_hogs_zero_or_one_refined"] HOGs changed #e58e73" >> $sankey_out
        echo "1-to-1 HOGs USC ["$number_of_hogs_one_to_one_USC"] HOGs USC" >> $sankey_out
        echo "Paralog HOGs USC ["$number_paralogs_USC"] HOGs USC" >> $sankey_out
        echo "0 or 1 gene HOGs USC ["$number_zero_or_one_USC"] HOGs USC" >> $sankey_out
        echo "1-to-1 HOGs USU ["$number_of_hogs_one_to_one_USU"] HOGs USU #8db600" >> $sankey_out
        echo "Paralog HOGs USU ["$number_paralogs_USU"] HOGs USU #8db600" >> $sankey_out
        echo "0 or 1 gene HOGs USU ["$number_zero_or_one_USU"] HOGs USU #8db600" >> $sankey_out
        echo ":Paralog HOGs USC #228B22" >> $sankey_out
        echo ":Paralog HOGs USU #228B22" >> $sankey_out
        echo ":1-to-1 HOGs USC #228B22" >> $sankey_out
        echo ":1-to-1 HOGs USU #228B22" >> $sankey_out
        echo ":0 or 1 gene HOGs USC #228B22" >> $sankey_out
        echo ":0 or 1 gene HOGs USU #228B22" >> $sankey_out
        echo ":1-to-1 HOGs CWS #ff0000" >> $sankey_out
        echo ":0 or 1 gene HOGs CWS #ff0000" >> $sankey_out
        echo ":Paralog HOGs CWS #ff0000" >> $sankey_out
    else  # complex == 4
        echo "Paralog HOGs ["$number_hogs_with_paralogs_refined"] Paralog HOGs changed #800080" >> $sankey_out
        number_paralogs_unchanged=$(($number_paralogs_USC+$number_paralogs_USU))
        echo "Paralog HOGs ["$number_paralogs_unchanged"] Paralog HOGs unchanged #800080" >> $sankey_out
        #echo "Paralog HOGs ["$number_paralogs_USU"] Paralog HOGs unchanged #800080" >> $sankey_out
        echo "1-to-1 HOGs ["$number_hogs_one_to_one_refined"] 1-to-1 HOGs changed" >> $sankey_out
        number_1_to_1_unchanged=$(($number_of_hogs_one_to_one_USC+$number_of_hogs_one_to_one_USU))
        echo "1-to-1 HOGs ["$number_1_to_1_unchanged"] 1-to-1 HOGs unchanged" >> $sankey_out
        #echo "1-to-1 HOGs ["$number_of_hogs_one_to_one_USU"] 1-to-1 HOGs unchanged" >> $sankey_out
        echo "0 or 1 gene HOGs["$number_hogs_zero_or_one_refined"] 0 or 1 gene HOGs changed" >> $sankey_out
        number_zero_or_one_unchanged=$(($number_zero_or_one_USC+$number_zero_or_one_USU))
        echo "0 or 1 gene HOGs["$number_zero_or_one_unchanged"] 0 or 1 gene HOGs unchanged" >> $sankey_out
        #echo "0 or 1 gene HOGs["$number_zero_or_one_USU"] 0 or 1 gene HOGs unchanged" >> $sankey_out
        echo "Paralog HOGs changed["$number_hogs_with_paralogs_refined"] HOGs changed #e58e73" >> $sankey_out
        echo "1-to-1 HOGs changed ["$number_hogs_one_to_one_refined"] HOGs changed #e58e73" >> $sankey_out
        echo "0 or 1 gene HOGs changed["$number_hogs_zero_or_one_refined"] HOGs changed #e58e73" >> $sankey_out
        echo "1-to-1 HOGs unchanged ["$number_1_to_1_unchanged"] HOGs unchanged" >> $sankey_out
        echo "Paralog HOGs unchanged ["$number_paralogs_unchanged"] HOGs unchanged" >> $sankey_out
        echo "0 or 1 gene HOGs unchanged ["$number_zero_or_one_unchanged"] HOGs unchanged" >> $sankey_out
        #echo "1-to-1 HOGs unchanged ["$number_of_hogs_one_to_one_USU"] HOGs unchanged #8db600" >> $sankey_out
        #echo "Paralog HOGs unchanged ["$number_paralogs_USU"] HOGs unchanged #8db600" >> $sankey_out
        #echo "0 or 1 gene HOGs unchanged ["$number_zero_or_one_USU"] HOGs unchanged #8db600" >> $sankey_out
        echo ":Paralog HOGs unchanged #228B22" >> $sankey_out
        echo ":1-to-1 HOGs unchanged #228B22" >> $sankey_out
        #echo ":1-to-1 HOGs USU #228B22" >> $sankey_out
        echo ":0 or 1 gene HOGs unchanged #228B22" >> $sankey_out
        #echo ":0 or 1 gene HOGs USU #228B22" >> $sankey_out
        echo ":1-to-1 HOGs changed #ff0000" >> $sankey_out
        echo ":0 or 1 gene HOGs changed #ff0000" >> $sankey_out
        echo ":Paralog HOGs changed #ff0000" >> $sankey_out
    fi

    if [ $complex == 2 ] || [ $complex == 3 ]; then
        echo "HOGs USC ["$number_hogs_USC"] HOGs unchanged #0F0" >> $sankey_out
        echo "HOGs USU ["$number_hogs_USU"] HOGs unchanged #0F0" >> $sankey_out
    fi

    echo ":HOGs changed #ff0000" >> $sankey_out
    echo ":Paralog HOGs #800080" >> $sankey_out
    echo ":HOGs unchanged #228B22" >> $sankey_out
    if [ $complex != 1 ]; then
        echo ":HOGs only 1 genome #000000" >> $sankey_out
    fi

    # figure settings
    echo "size w 1200" >> $sankey_out
    echo "  h 1200" >> $sankey_out
    echo "layout order automatic" >> $sankey_out
    echo "  justifyends Y" >> $sankey_out

    if [ $complex == 2 ] || [ $complex == 3 ]; then
        echo "
        1-to-1 HOGs are HOGs that contained precisely one gene per genome. 
        0 or 1 gene HOGs are HOGs missing an ortholog in at least one genome and none of the genomes contributed more than one gene. 
        Paralog HOGs are the HOGs where at least one genome contributed more than one gene. Unchanged synteny confirmed (USC) are the HOGs that remained the same after OrthoRefine and the groupings had synteny support. 
        Unchanged synteny unconfirmed (USU) are the HOGs that remained the same after OrthoRefine, but the groups did not have synteny support. 
        Changed with synteny (CWS) are the HOGs where a refinement could be made as supported by synteny. HOGs where only one genome contributed several genes, and no other genomes were present, can’t be analyzed with OrthoRefine."
    fi
fi

echo "Total HOGs: $total_hogs
Total 1-to-1 HOGs: $one_to_one_orthologs
Total 0-or-1 HOGs: $zero_or_one_hogs
Total Paralog HOGs: $hogs_with_paralogs
Confirmed 1-to-1 HOGs: $number_of_hogs_one_to_one_USC
Unconfirmed 1-to-1 HOGs: $number_of_hogs_one_to_one_USU
Modified 1-to-1 HOGs: $number_hogs_one_to_one_refined
Confirmed 0-or-1 HOGs: $number_zero_or_one_USC
Unconfirmed 0-or-1 HOGs: $number_zero_or_one_USU
Modified 0-or-1 HOGs: $number_hogs_zero_or_one_refined
Confirmed Paralog HOGs: $number_paralogs_USC
Unconfirmed Paraolg HOGs: $(( $number_paralogs_USU+$number_hogs_only_1_genome_no_refine ))
Modified Paralog HOGs: $number_hogs_with_paralogs_refined
Total Confirmed HOGs: $(( $number_of_hogs_one_to_one_USC+$number_zero_or_one_USC+$number_paralogs_USC ))
Total Unconfirmed HOGs: $(( $number_of_hogs_one_to_one_USU+$number_zero_or_one_USU+$number_paralogs_USU+$number_hogs_only_1_genome_no_refine ))
Total Modified HOGs: $(( $number_hogs_one_to_one_refined+$number_hogs_zero_or_one_refined+$number_hogs_with_paralogs_refined ))"

echo "Total HOGs: $total_hogs
Total 1-to-1 HOGs: $one_to_one_orthologs
Total 0-or-1 HOGs: $zero_or_one_hogs
Total Paralog HOGs: $hogs_with_paralogs
Confirmed 1-to-1 HOGs: $number_of_hogs_one_to_one_USC
Unconfirmed 1-to-1 HOGs: $number_of_hogs_one_to_one_USU
Modified 1-to-1 HOGs: $number_hogs_one_to_one_refined
Confirmed 0-or-1 HOGs: $number_zero_or_one_USC
Unconfirmed 0-or-1 HOGs: $number_zero_or_one_USU
Modified 0-or-1 HOGs: $number_hogs_zero_or_one_refined
Confirmed Paralog HOGs: $number_paralogs_USC
Unconfirmed Paraolg HOGs: $(( $number_paralogs_USU+$number_hogs_only_1_genome_no_refine ))
Modified Paralog HOGs: $number_hogs_with_paralogs_refined
Total Confirmed HOGs: $(( $number_of_hogs_one_to_one_USC+$number_zero_or_one_USC+$number_paralogs_USC ))
Total Unconfirmed HOGs: $(( $number_of_hogs_one_to_one_USU+$number_zero_or_one_USU+$number_paralogs_USU+$number_hogs_only_1_genome_no_refine ))
Total Modified HOGs: $(( $number_hogs_one_to_one_refined+$number_hogs_zero_or_one_refined+$number_hogs_with_paralogs_refined ))"> OrthoRefine_summary.txt