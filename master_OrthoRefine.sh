#!/usr/bin/env bash

finder_cmd=""
refine_cmd=""
OF_file="N0.tsv"

blast_stop=false
othrogroups_stop=false
sequence_stop=false
alignment_stop=false
tree_stop=false
orthogroups_start=false
tree_start=false

while true; do
    thing=$1
    if [ -z $thing ];then #if no more arguments can be loaded
        shift
        break
    elif [ $thing == --OrthoFinder ];then
        finder_cmd=$2$finder_cmd
        shift 2
    elif [ $thing == --OrthoRefine ];then
        refine_cmd=$2$refine_cmd
        shift 2
    elif [ $thing == -op ] && [ $blast_stop == false ];then
        blast_stop=true
        finder_cmd=$finder_cmd" -op"
        shift
    elif [ $thing == -og ] && [ $orthogroups_stop == false ];then
        echo "test"
        othrogroups_stop=true
        finder_cmd=$finder_cmd" -og"
        shift
    elif [ $thing == -os ] && [ $sequence_stop == false ];then
        sequence_stop=true
        finder_cmd=$finder_cmd" -os"
        shift
    elif [ $thing == -oa ] && [ $alignment_stop == false ];then
        alignment_stop=true
        finder_cmd=$finder_cmd" -oa"
        shift
    elif [ $thing == -ot ] && [ $tree_stop == false ];then
        tree_stop=true
        finder_cmd=$finder_cmd" -ot"
        shift
    elif [ $thing == -fg ] && [ $orthogroups_start == false ];then
        orthogroups_start=true
        finder_cmd=$finder_cmd" -fg"
        shift
    elif [ $thing == -ft ] && [ $tree_start == false ];then
        tree_start=true
        finder_cmd=$finder_cmd" -ft"
        shift
    elif [ $thing == -f ];then
        fasta="$2"
        finder_cmd=$finder_cmd" -f $fasta"
        shift 2
    elif [ $thing == -b ];then
        precompute_blast="$2"
        finder_cmd=$finder_cmd" -b $precompute_blast"
        shift 2
    elif [ $thing == -t ];then
        search_threads="$2"
        finder_cmd=$finder_cmd" -t $search_threads"
        shift 2
    elif [ $thing == -a ];then
        analysis_threads="$2"
        finder_cmd=$finder_cmd" -a $analysis_threads"
        shift 2
    elif [ $thing == -d ];then
        input_DNA=true
        finder_cmd=$finder_cmd" -d"
        shift
    elif [ $thing == -M ];then
        method_gene_tree_inf="$2".
        finder_cmd=$finder_cmd" -M $method_gene_tree_inf"
        shift 2
    elif [ $thing == -S ];then
        search_prog="$2"
        finder_cmd=$finder_cmd" -S $search_prog"
        shift 2
    elif [ $thing == -A ];then
        MSA_prog="$2"
        finder_cmd=$finder_cmd" -A $MSA_prog"
        shift 2
    elif [ $thing == -T ];then
        method_tree_inf="$2"
        finder_cmd=$finder_cmd" -T $method_tree_inf"
        shift 2
    elif [ $thing == -s ];then
        user_root_tree="$2"
        finder_cmd=$finder_cmd" -s $user_root_tree"
        shift 2
    elif [ $thing == -I ];then
        MCL_inflat="$2"
        finder_cmd=$finder_cmd" -I $MCL_inflat"
        shift 2
    elif [ $thing == -x ];then
        XML_form="$2"
        finder_cmd=$finder_cmd" -x $XML_form"
        shift 2
    elif [ $thing == -p ]; then
        pickle_file="$2"
        finder_cmd=$finder_cmd" -p $pickle_file"
        shift 2
    elif [ $thing == -1 ];then
        one_way_search=true
        finder_cmd=$finder_cmd" -1"
        shift
    elif [ $thing == -X ];then
        spec_name=true
        finder_cmd=$finder_cmd" -X"
        shift
    elif [ $thing == -y ];then
        split_paralog=true
        finder_cmd=$finder_cmd" -y"
        shift
    elif [ $thing == -z ];then
        no_trim_MSA=true
        finder_cmd=$finder_cmd" -z"
        shift
    elif [ $thing == -n ];then
        append_name="$2"
        finder_cmd=$finder_cmd" -n $append_name"
        shift 2
    elif [ $thing == -o ];then
        res_dir="$2" # don't remove this line
        finder_cmd=$finder_cmd" -o $res_dir"
        shift 2
    elif [ $thing == --input ];then
        input="$2" # input file with REFSEQ GCF for OrthoRefine
        refine_cmd=$refine_cmd" --input $input"
        shift 2
    elif [ $thing == --path ];then
        path="$2"
        refine_cmd=$refine_cmd" --path $path"
        shift 2
    elif [ $thing == --window_size ];then
        window_size="$2"
        refine_cmd=$refine_cmd" --window_size $window_size"
        shift 2
    elif [ $thing == --synteny_ratio ];then
        synteny_ratio=$2
        refine_cmd=$refine_cmd" --synteny_ratio $synteny_ratio"
        shift 2
    elif [ $thing == --diag ];then
        diag="$2"
        refine_cmd=$refine_cmd" --diag $diag"
        shift 2
    elif [ $thing == --OF_file ];then
        OF_file="$2"
        refine_cmd=$refine_cmd" --OF_file $OF_file"
        shift 2
    elif [ $thing == --outfile ];then
        outfile="$2"
        refine_cmd=$refine_cmd" --outfile $outfile"
        shift 2
    elif [ $thing == --print_all ];then
        print_all="$2"
        refine_cmd=$refine_cmd" --print_all $print_all"
        shift 2
    elif [ $thing == --run_all_orthofinder ];then
        run_all_orthofinder="$2"
        refine_cmd=$refine_cmd" --run_all_orthofinder $run_all_orthofinder"
        shift 2
    elif [ $thing == --benchmark ];then
        benchmark="$2"
        refine_cmd=$refine_cmd" --benchmark $benchmark"
        shift 2
    elif [ $thing == --run_single_HOG ];then
        run_single_Hog="$2"
        refine_cmd=$refine_cmd" --run_single_HOG $run_single_Hog"
        shift 2
    elif [ $thing == --paralogs_print ];then
        paralogs_print="$2"
        refine_cmd=$refine_cmd" --paralogs print $paralogs_print"
        shift 2
    elif [ $thing == --prod_acc ];then
        prod_acc="$2"
        refine_cmd=$refine_cmd" --prod_acc $prod_acc"
    else
    echo $@
        echo "Error: bad argument provided " $thing
        exit
    fi
done

./download_ft_fafiles.sh $input

$finder_cmd
#echo $finder_cmd

if [ ! -z "$res_dir" ];then # if res_dir was set by user
    cp "$(\ls -1dt $res_dir/OrthoFinder/*/Phylogenetic_Hierarchical_Orthogroups/$OF_file | head -n 1)" ./
else
   cp "$(\ls -1dt ./OrthoFinder/*/Phylogenetic_Hierarchical_Orthogroups/$OF_file | head -n 1)" ./
fi

dos2unix $OF_file

$refine_cmd
echo $refine_cmd

