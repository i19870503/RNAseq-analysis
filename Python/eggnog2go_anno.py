import os 
import re
import pandas as pd
import string
import itertools
import numpy as np
import sys
import argparse
from collections import OrderedDict

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create GO annotation and enrichment file')
    parser.add_argument('-i',type=str,dest='infile',required=True,help="Input file")
    parser.add_argument('-o',type=str,dest='out',required=True,help="Ouput file")
    parser.add_argument('-db',type=str,dest='db',required=True,help="GO Database file")
    args = parser.parse_args()
    print (args)

def sort_uniq(sequence):
    return (x[0] for x in itertools.groupby(sorted(sequence)))

path = "/home/zluna/Work/GO"
fout = open(args.out+"_anno.xls", 'w')
print("Gene_id", "GO_annotation", sep = '\t', file = fout)
go_db = pd.read_table(os.path.join(path, args.db), header = None)
eggout = pd.read_table(os.path.join(path, args.infile), header = None)
#pd.DataFrame.head(eggout)
#eggout.head(100)
dict = OrderedDict()
first_flag = 1 
a = list(go_db[0])
for i in range(len(eggout)):
    gene_id = eggout[0][i]
    go_id = eggout[5][i]
    if pd.isnull(eggout[5][i]):
        go_id = ''
    #print(gene_id, kegg_id, type(kegg_id), sep ='\t')
    go_id = go_id.split(',')
    if len(go_id) == 0:
        continue
    go_term = '; '.join(list(go_db[go_db[2].isin(go_id)][0]))
    #print(gene_id, go_id, go_term, sep ='\t')
    go_sum = []
    sel_go_table = go_db[go_db[2].isin(go_id)]
    for j in range(len(sel_go_table)):
        go_sum.append(''.join(( list(sel_go_table[2])[j], "~", list(sel_go_table[0])[j])))
    print(gene_id, str(go_sum).strip('[]').replace(']','').replace("'","").replace(", ","; "), sep = '\t', file = fout)
    

    a = list(go_db[2])
### Use dictionary 
    for k in range(len(a)):
        if str(go_sum).find(a[k]) != -1 :
            if a[k] not in dict.keys():
### The value must be list type, if just give the 'gene_id' as the value of key, it can not use 'append' method to add the new 'gene_id' to the existing key.
                dict[a[k]] = []
                dict[a[k]].append(gene_id)
            else:                                
                dict[a[k]].append(gene_id)
                
                #dict[a[j]] = [dict[a[j]], gene_id]
fout.close()

fout2 = open(args.out+"_enrich.xls", 'w')        
print('GOID', 'Term', 'Genes', 'Gene_count', sep = '\t', file = fout2)
for key,values in  dict.items():
    print(key, list(go_db[go_db[2] == key][0]), str(values).strip('[]').replace(']','').replace("'",""), len(values), sep ='\t', file = fout2)
fout2.cloes()
