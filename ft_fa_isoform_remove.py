# Script to remove isoforms from feature table file based on either first occurance (+) strand OR last occurance (-) strand
# Script modifies both feature table and fasta file. The fasta file has the identified duplicates, from the feature table, removed.
# Input is GCF file prefix. E.g. GCF_000001405.40
# Script by J. Ludwig. 2023. 

import os
import fnmatch
import sys
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def get_file_name(file_prefix):
    for file in os.listdir('.'):
        if fnmatch.fnmatch(file, file_prefix):
            return file
    return None

def fasta_to_df(file_name):
    iden = []
    seqs = []
    desc = []
    for entry in SeqIO.parse(file_name, "fasta"):
        iden.append(entry.id)
        seqs.append(entry.seq)
        desc.append(entry.description)
    return pd.DataFrame({'id': iden, 'seq': seqs, 'desc': desc})

def df_to_fasta(df, filename, line_length=80):
    with open(filename, 'w') as f:
        for index, row in df.iterrows():
            seq = Seq(row['seq'])
            record = SeqRecord(seq, id=row['id'])
            record.description = row['desc']
            f.write(custom_seq_record_to_string(record, line_length) + "\n")

def custom_seq_record_to_string(seq_record, line_length=80):
    # Convert the SeqRecord to a string in FASTA format
    lines = []
    lines.append(">" + seq_record.description)
    for i in range(0, len(seq_record.seq), line_length):
        lines.append(str(seq_record.seq[i:i+line_length]))
    return "\n".join(lines)


file_prefix = sys.argv[1]
# isoform_pref = sys.argv[2] # not implemented yet. Will be first occurance or longest sequence

ft_file_name = get_file_name(file_prefix + '*_feature_table.txt')
ft_df = pd.read_csv(ft_file_name, sep='\t', dtype = str)
#only keep rows that have "CDS" in the #feature column
ft_df = ft_df[ft_df['# feature'].str.contains("CDS")]
#only keep rows that do not have "without" in the class column
ft_df = ft_df[~ft_df['class'].str.contains("without")]
# if on + strand, keep the first duplicate, if on - strand, keep the last duplicate
ft_df['duplicate_for'] = ft_df.loc[ft_df['strand'] == '+'].duplicated(subset=['symbol'], keep='first')
ft_df['duplicate_rev'] = ft_df.loc[ft_df['strand'] == '-'].duplicated(subset=['symbol'], keep='last')
ft_df.drop(ft_df.loc[ft_df['duplicate_for'] == True].index, inplace=True)
ft_df.drop(ft_df.loc[ft_df['duplicate_rev'] == True].index, inplace=True)
# have to redrop duplicates as the same symbol can be on both strands (alternative scaffolds)
ft_df.drop_duplicates(subset=['symbol'], keep='first', inplace=True)

print(ft_df.iloc[:, 0:20].to_csv(ft_file_name, sep = '\t', index = False))

fa_file_name = get_file_name(file_prefix + '*_protein.faa')
#read fasta file into dataframe
fa_df = fasta_to_df(fa_file_name)
#only keep rows that are in ft_df
fa_df = fa_df[fa_df['id'].isin(ft_df['product_accession'])]
#print(fa_df.to_csv(fa_test, sep = '\t', index = False))
df_to_fasta(fa_df, fa_file_name)

