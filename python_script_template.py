#!/usr/bin/python
# -*- coding: utf-8 -*-

## Function: Create files_summary_virus.txt for building virus database

## Verison log:
# 2021-08-31 Created

## Usage: python create_files_summary_virus.py \
#         --inputfile APG_NCBIVirusDB_merged_yyyy-mm-dd.csv \
#         --genomePath ./genome \
#         --outpath ./output 
# 
# input:(1) APG_NCBIVirusDB_merged_2021-08-25.csv
# output:(2) files_summary_virus.txt

#python package
import os
import time
import argparse
import pandas as pd
import numpy as np
import csv
from pandas._libs.missing import NA

## time start record
d1 = time.time()

print('____________________________')

########################
# commend line arguments
########################

# Usage
message_prog = 'python create_files_summary_virus.py'
message_description = 'Create files_summary_virus.txt for building virus database'
# example
message_epilog = """
python3 create_files_summary_virus.py \
--inputfile APG_NCBIVirusDB_merged_yyyy-mm-dd.csv \
--genomePath ./genome \
--outpath ./output 
"""

## set argparse
parser = argparse.ArgumentParser(
    prog=message_prog, # usage
    description=message_description, # function description
    epilog=message_epilog)


parser.add_argument('-i','--inputfile', type=str, help='APG_NCBIVirusDB_merged_yyyy-mm-dd.csv, which is created by APG_NCBIVirusDB.py')
parser.add_argument('-gp','--genomePath', type=str, help='the path containing virus genome')
parser.add_argument('-o','--outpath', default='./', type=str, help='Output text file path')



# get parse_args object
args = parser.parse_args()

## put commend line into variable
input_file = args.inputfile
genome_path = args.genomePath
path_output = args.outpath

##########################
# Path check
##########################

if genome_path is None:
    print('ERROR: genome path empty')
    parser.print_help()
    quit()


# error check input file exist or not
if os.path.exists(input_file):
    print('The file: ' + input_file + ' exist')
    print('____________________________')
else:
    print('ERROR:' + input_file + ' file not exist')
    quit()

if os.path.exists(genome_path):
    print('The directory: ' + genome_path + ' exist')
    print('____________________________')
else:
    print('ERROR:' + genome_path + ' file not exist')
    quit()

if genome_path[-1] == '/':
    print('Input genome path correct')
else:
    genome_path = genome_path + '/'

##########################
# output dir 
##########################
# 
if os.path.exists(path_output):
    pass
else:
    os.system('mkdir ' + path_output )
    print('Creat New Directory:' + path_output)

if path_output[-1] == '/':
    print('Output path is OK.')
else:
    path_output = path_output + '/'


##########################
# Star to process 
##########################

# extract species names
NCBI_table = pd.read_csv(input_file, sep = ',')
print('The number of species is', NCBI_table.shape[0])
NCBI_table['Species'].to_csv(path_output + 'virus_list.txt', sep = '\t', index=False, header=False)
print('____________________________')

# get taxid by taxonkit
os.system('rm -f ' + path_output + 'names2taxid.txt')
os.system('cat ' + path_output + 'virus_list.txt | taxonkit name2taxid | taxonkit lineage -i 2 -r -L >> ' + path_output + 'names2taxid.txt')
print('____________________________')


# extract the species without taxid
os.system('rm -f ' + path_output + 'sp_without_taxid.txt')
os.system("awk -F '\t' '$2 !~ /[0-9]/' " + path_output + "names2taxid.txt | uniq >> " + path_output + "sp_without_taxid.txt")
os.system("sed -i 's/\t//g' " + path_output + "sp_without_taxid.txt") # remove \t
os.system("echo 'The number of species without taxid is '$(cat " + path_output + "sp_without_taxid.txt | wc -l)")
print('____________________________')

# get taxid artifiically by greping names.dmp
if os.path.exists(path_output + 'sp_still_miss.csv') and os.path.exists(path_output + 'sp_without_taxid_grep.txt'):

    sp_manually_added = pd.read_csv(path_output + 'sp_still_miss.csv', sep = ',', header=None)
    print('The number of manually-added species is ' +  str(len(sp_manually_added[0])))

    print('____________________________')
    sp_without_taxid_grep = pd.read_csv(path_output + 'sp_without_taxid_grep.txt', sep = '\t', header=None)
    print('Add sp_still_miss.csv to sp_without_taxid_grep.txt.')

    #sp_without_taxid_grep = sp_without_taxid_grep.iloc[:,[0,2]] # select column
    #sp_without_taxid_grep = sp_without_taxid_grep.set_axis(['0', '1'], axis=1) # rename col name
    #sp_manually_added = sp_manually_added.set_axis(['0', '1'], axis=1) # rename col name
    sp_merged = pd.concat([sp_without_taxid_grep,sp_manually_added], ignore_index=True) # merge table
    sp_merged.to_csv(path_output + 'sp_without_taxid_grep.txt', sep = '\t', index=False, header=False) # write merged table
    

else:
    print('Grepping ...')
    os.system("IFS=$'\n'; rm " + path_output + "sp_without_taxid_grep.txt; for i in $(cat " + path_output + "sp_without_taxid.txt); do echo $i ; cat ~/.taxonkit/names.dmp | grep -e $i | awk -F '|' '{print $1, $2}' >> " + path_output + "sp_without_taxid_grep.txt; done")
    os.system('if [ -s "' + path_output + 'sp_without_taxid_grep.txt" ];then echo "Some species were catched.";else echo "Nothing is catched.";exit;fi')

    # sp with still no taxid
    sp_without_taxid = pd.read_csv(path_output + 'sp_without_taxid.txt', sep = '\t', header=None)
    sp_without_taxid_grep = pd.read_csv(path_output + 'sp_without_taxid_grep.txt', sep = '\t', header=None)
    sp_still_no_taxid = set(sp_without_taxid[0]).difference(sp_without_taxid_grep[2])
    
    if len(sp_still_no_taxid)>0:
        print('The taxid of following species are still not detected: \n' + str(sp_still_no_taxid))
        print('Total: ' + str(len(sp_still_no_taxid)))
        print('Please add the taxid of these species manually.')
        print('Create the file "sp_still_miss.csv", the content: \n' + '1111,taxon1 \n2222,taxon2 \n' + '. \n' + '. \n' + '. \n')
        print('After finishing the file, run the script again.')
        quit()
    else:
        print('Great! There is no species without taxid.')

# exactly extract
print('Exactly extract')
sp_without_taxid = pd.read_csv(path_output + 'sp_without_taxid.txt', sep = '\t', header=None)
sp_without_taxid_grep = pd.read_csv(path_output + 'sp_without_taxid_grep.txt', sep = '\t', header=None)
sp_without_taxid_grep_exac = sp_without_taxid_grep[sp_without_taxid_grep[1].isin(sp_without_taxid[0])]
sp_without_taxid_grep_exac = sp_without_taxid_grep_exac.drop_duplicates()
sp_without_taxid_grep_exac.to_csv(path_output + 'sp_without_taxid_grep_exac.txt', sep = '\t', index=False, header=False)

# make files_summary.txt
print('Making files_summary.txt')
sp_taxid = pd.read_csv(path_output + 'names2taxid.txt', sep='\t', header=None, dtype=str)
sp_with_taxid = sp_taxid.iloc[:,[0,1]].dropna()
sp_without_taxid_grep_exac.iloc[:, [0,1]] = sp_without_taxid_grep_exac.iloc[:, [1,0]]
sp_all_taxid_id = sp_with_taxid.append(sp_without_taxid_grep_exac, ignore_index=True)

# NCBI_table.replace(sp_all_taxid_id[0].tolist(), sp_all_taxid_id[1].tolist())

d = {'file_name': NCBI_table['Accession'],
      'source': ['NCBIVirus']*NCBI_table.shape[0],
      'Assembly_accession': NCBI_table['Accession'],
      'type': ['Virus']*NCBI_table.shape[0],
      'Species_taxid': NCBI_table['Species'],
      'taxid': NCBI_table['Species'],
      'species_name': NCBI_table['Species'],
      'Geo_Location': NCBI_table['Geo_Location'],
      'Segment': NCBI_table['Segment'],
      'Sequence_Type': NCBI_table['Sequence_Type']}

df = pd.DataFrame(data=d)
df['Species_taxid'] = df['Species_taxid'].replace(sp_all_taxid_id[0].tolist(), sp_all_taxid_id[1].tolist())
df['taxid'] = df['taxid'].replace(sp_all_taxid_id[0].tolist(), sp_all_taxid_id[1].tolist())
df['file_name'] = [x + '.fna.gz' for x in df['file_name']]
df.to_csv(path_output + 'files_summary.txt', sep='\t', index=False, header=True)


# get the segment and merge the .fa
print('Get the segment and merge the .fa')
df_Refseq = df[df['Sequence_Type'].isin(['RefSeq'])]
df_Refseq_seg = df_Refseq[~df_Refseq['Segment'].isin([np.NAN])]
segment_sp = df_Refseq_seg['species_name']
segment_fa = [x + '.fa' for x in df_Refseq_seg['Assembly_accession']]

d = {'species_name': segment_sp,
     'file_name': segment_fa}

d_seg = pd.DataFrame(data=d)
d_seg.to_csv(path_output + 'sp_segment.txt', sep='\t', index=False, header=True)

print('____________________________')
print('Create segment genome path')
os.system('mkdir ' + path_output + 'seg_merge_genome; rm ' + path_output + 'seg_merge_genome/*')


df_new = df
sp_n = len(d_seg['species_name'].drop_duplicates())
k=0
for i in d_seg['species_name'].drop_duplicates():
    k += 1
    print(i)
    for j in d_seg[d_seg['species_name'].isin([i])]['file_name']:
        print(j)
        print(k)
        if os.path.exists(genome_path + j ):
            os.system('cat ' + genome_path + j + ' >> ' + path_output + 'seg_merge_genome/APG_seg_' + str(k) + '.fa')
            n = len(df_new['file_name'][df_new['species_name'].isin([i])])
            df_new['file_name'][df_new['species_name'].isin([i])] = ['APG_seg_' + str(k) + '.fna.gz']*n
            print(str(k) + '/' + str(sp_n))
        else:
            print(genome_path + j + ' not found')
            quit()

if os.path.exists(path_output + 'seg_merge_genome'):
    os.system('gzip ' + path_output + 'seg_merge_genome/*.fa ; rename .fa.gz .fna.gz ' + path_output + 'seg_merge_genome/*.fa.gz')
    os.system('gzip ' + genome_path + '*.fa ; rename .fa.gz .fna.gz ' + genome_path + '*.fa.gz; cp ' + path_output + 'seg_merge_genome/* ' + genome_path)
    print('____________________________')
    print('Merged segment genome has been copied to genome path')
else:
    print('The directory: ' + path_output + 'seg_merge_genome, is not found')
    quit()


print('____________________________')
print('Making final result: files_summary_virus.txt')
df_new['Species_taxid'] = df_new['Species_taxid'].replace(sp_all_taxid_id[0].tolist(), sp_all_taxid_id[1].tolist())
df_new['taxid'] = df_new['taxid'].replace(sp_all_taxid_id[0].tolist(), sp_all_taxid_id[1].tolist())
df_new.to_csv(path_output + 'files_summary_virus.txt', sep='\t', index=False, header=True)

