#!/bin/bash
IFS=$'\n | \t' read -d '' -a DOWNLOAD < $1                                      # Read GCF accession from this file. One sample per line. Store in array "DOWNLOAD".

length=$(( ${#DOWNLOAD[@]} - 1 ))                                               # Number of samples from file.

cl_flag=0                                                                       # flag to check if second column contains "c" "C" "l" "L" for circular or linear genome
for j in `seq 0 $length`; do
  if [ "${DOWNLOAD[j]}" == c ] || [ "${DOWNLOAD[j]}" == C ] || [ "${DOWNLOAD[j]}" == l ] || [ "${DOWNLOAD[j]}" == L ]; then
    cl_flag=1
    break
  fi
done

abe_flag=0                                                                      # archaea, bacteria, eukaryote flag
for j in `seq 0 $length`; do
  if [ "${DOWNLOAD[j]}" == a ] || [ "${DOWNLOAD[j]}" == A ] || [ "${DOWNLOAD[j]}" == b ] || [ "${DOWNLOAD[j]}" == B ] || [ "${DOWNLOAD[j]}" == e ] || [ "${DOWNLOAD[j]}" == E ]; then
    abe_flag=1
    break
  fi
done


if [ $cl_flag == 1 ]; then
  if ! (($length%2)) && ! (($length+1%2)); then                                                      # checks if any second column has "c" "C" "l" "L" and if any second column is missing data
	  echo "ERROR: Odd number of inputs from file. Check each column has an input"
    return 1 2>/dev/null
    { exit 1; }
  fi
  if [ $abe_flag == 0 ]; then
    for j in `seq 1 2 $length`; do                                                # checks if all second columns contain a letter and that letter is "c" "C" "l" "L"
      if ! [ "${DOWNLOAD[j]}" == c ] && ! [ "${DOWNLOAD[j]}" == C ] && ! [ "${DOWNLOAD[j]}" == l ] && ! [ "${DOWNLOAD[j]}" == L ]; then
        echo "ERROR: Second column contains not "c" or "C" or "l" or "L""
        return 1 2>/dev/null
        { exit 1; }
      fi
    done
  else
    for j in `seq 1 3 $length`; do                                                # checks if all second columns contain a letter and that letter is "c" "C" "l" "L"
      if ! [ "${DOWNLOAD[j]}" == c ] && ! [ "${DOWNLOAD[j]}" == C ] && ! [ "${DOWNLOAD[j]}" == l ] && ! [ "${DOWNLOAD[j]}" == L ]; then
        echo "ERROR: Second column contains not "c" or "C" or "l" or "L""
        return 1 2>/dev/null
        { exit 1; }
      fi
    done
    for j in `seq 2 3 $length`; do
      if ! [ "${DOWNLOAD[j]}" == a ] && ! [ "${DOWNLOAD[j]}" == A ] && ! [ "${DOWNLOAD[j]}" == b ] && ! [ "${DOWNLOAD[j]}" == B ] && ! [ "${DOWNLOAD[j]}" == e ] && ! [ "${DOWNLOAD[j]}" == E ]; then
      echo "ERROR: Third column contains not "a" or "A" or "b" or "B" or "e" or "E""
      return 1 2>/dev/null
      { exit 1; }
      fi
    done
  fi
fi



motd_flag=0                                                                     # Flag so MOTD from NCBI only printed once.
j=0
while [ $j -le $length ]; do                                                    # Loop through elements of "DOWNLOAD".
  check=${DOWNLOAD[j]}                                                          # # Check if both protein and feature_table file are already downloaded and then if so, skip downloading them again
  if [ -f $check\_*_feature_table.txt ] && [ -f $check\_*_protein.faa ]; then
    ((j+=1))                                                                      # increment every time loop runs
    if [ $cl_flag == 1 ]; then                                                    # increment again if there are 2 columns present as we need to skip the second column data
      ((j+=1))
    fi
    if [ $abe_flag == 1 ]; then
      ((j+=1))
    fi
    continue
  fi

  GCX=${DOWNLOAD[j]:0:3}                                                        # GCF or GCA
  add=                                                                          # Variable to store every 3 numbers after "GCF_" and before ".#" of accession.
  length2=${#DOWNLOAD[j]}                                                       # Number of characters in array element.
  for i in `seq 4 $length2`; do                                                 # Start at first number. Skip "GCF_".
    if [ "${DOWNLOAD[j]:$i:1}" == "." ]; then                                   # Stop storing characters (numbers) at first "."
      add=$add${DOWNLOAD[j]}                                                    # Append full GCF accession from user input; used for when mutiple versions are present on FTP; only download desired version. GCF_1234.1 v GCF_1234.2
      break
    fi
    if [ "$(( $i % 3))" == "0" ]; then                                          # Store every 3 numbers followed by "/".
      add=$add${DOWNLOAD[j]: $(( $i - 2 )) :3}/
    fi
  done

  add=rsync://ftp.ncbi.nlm.nih.gov/genomes/all/$GCX/$add*/                      # See end.
  if [ "$motd_flag" == "0" ]; then                                              # Only print NCBI MOTD once.
    motd_flag=1
    rsync -r --include "*_feature_table.txt.gz" --include "*_protein.faa.gz" --exclude="*" $add ./
  else
    rsync -r --no-motd --include "*_feature_table.txt.gz" --include "*_protein.faa.gz" --exclude="*" $add ./
  fi

  #old rsync code using exclude isntead of include
  #add=rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/$add*/                       # See below.
  #if [ "$motd_flag" == "0" ]; then                                              # Only print NCBI MOTD once.
  #  motd_flag=1
  #  rsync -r --exclude=*gpff* --exclude=md5checksums.txt --exclude=*rna*\
  #  --exclude=*annotation* --exclude=*genomic* --exclude=*assembly* \
  #  --exclude=*count* --exclude=*transl* --exclude=READ* --exclude=*wgsmaster* $add ./
  #else
  #  rsync -r --no-motd --exclude=*gpff* --exclude=md5checksums.txt --exclude=*rna*\
  #  --exclude=*annotation* --exclude=*genomic* --exclude=*assembly* \
  #  --exclude=*count* --exclude=*transl* --exclude=READ* --exclude=*wgsmaster* $add ./
  #fi

  check=${DOWNLOAD[j]}                                                          # Check if both protein and feature_table file were downloaded.
  if [ ! -f $check\_*_feature_table.txt.gz ]; then
    echo "Error:" $check\_*_feature_table.txt.gz "not downloaded"
  fi
  if [ ! -f $check\_*_protein.faa.gz ]; then
    echo "Error:" $check\_*_protein.faa.gz "not downloaded"
  fi

  gunzip $check\_*_feature_table.txt.gz                                         # unzip the downloaded files
  gunzip $check\_*_protein.faa.gz



  ((j+=1))                                                                      # increment every time loop runs
  if [ $cl_flag == 1 ]; then                                                    # increment again if there are 2 columns present as we need to skip the second column data
    ((j+=1))
  fi
  if [ $abe_flag == 1 ]; then
    ((j+=1))
  fi

done

# If the sample REFSEQ accession is "GCF_000009085.1"
# The two files we need and their web address are:
#   GCF_000009085.1_ASM908v1_protein.faa.gz
#     ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/085/GCF_000009085.1_ASM908v1/GCF_000009085.1_ASM908v1_protein.faa.gz
#   GCF_000009085.1_ASM908v1_feature_table.txt.gz
#     ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/085/GCF_000009085.1_ASM908v1/GCF_000009085.1_ASM908v1_feature_table.txt.gz
# Note how the accession is divided into 3 numbers: /GCF/000/009/085/
# The variable "add" will be "000/009/085/" at the end of the second for loop.
# The variable "add" is then combined with the static first part of the web address, "rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/".
# (rsync replaces ftp at the start because we will use rsync.)
# The variable "add" has the wildcard(*) and "/" added so the directory is not copied; just the two files we want inside it will be.
# (wildcard because the user is not required to include assembly version; in this example this is "ASM908v1".)
# The rsync command is issued:
# rsync -r --exclude=*gpff* --exclude=md5checksums.txt --exclude=*annotation* --exclude=*genomic* --exclude=*assembly* --exclude=*count*
#  --exclude=*transl* --exclude=READ* rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/009/085/*/
# (--exclude excludes the other files, based on their name, that we don't want to download)
