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
    parser = argparse.ArgumentParser(description='Create kegg annotation and enrichment file')
    parser.add_argument('-i',type=str,dest='infile',required=True,help="Input file")
    parser.add_argument('-o',type=str,dest='out',required=True,help="Ouput file")
    parser.add_argument('-db',type=str,dest='db',required=True,help="Database file")
    args = parser.parse_args()
    #print (args)

def sort_uniq(sequence):
    return (x[0] for x in itertools.groupby(sorted(sequence)))

path = "/home/zluna/Work/kegg"
fout = open(args.out, 'w')

print('Pathway', 'Genes', 'KOID', 'Gene_count', sep = '\t', file = fout)


df = pd.read_table(os.path.join(path, args.db))

eggout = pd.read_table(os.path.join(path, args.infile), header = None)
#pd.DataFrame.head(eggout)
#eggout.head()

dict = OrderedDict()
dict_koid = OrderedDict()
first_flag = 1 
for i in range(len(eggout)):
    gene_id = eggout[0][i]
    kegg_id = eggout[6][i]
    if pd.isnull(eggout[6][i]):
        kegg_id = ''
    #print(gene_id, kegg_id, type(kegg_id), sep ='\t')
    kegg_id = kegg_id.split(',')
    #print(kegg_id)
    kegg_path = '; '.join(list(df[df.KOID.isin(kegg_id)].Pathway))
    #print(gene_id, kegg_id, kegg_path, sep ='\t')
    
    a = list(filter(lambda x: re.search(r'PATH:', x), list(np.unique(list(df.Pathway)))))

    
### Use dictionary 
    for j in range(len(a)):
        if kegg_path.find(a[j]) != -1 :
            if a[j] not in dict.keys():
### The value must be list type, if just give the 'gene_id' as the value of key, it can not use 'append' method to add the new 'gene_id' to the existing key.
                dict[a[j]] = []
                dict_koid[a[j]] = []
                dict[a[j]].append(gene_id)
                dict_koid[a[j]].append(kegg_id)
            else:                                
                dict[a[j]].append(gene_id)
                dict_koid[a[j]].append(kegg_id)
                
                #dict[a[j]] = [dict[a[j]], gene_id]
                
for key,values in  dict.items():
    print(key, str(values).strip('[]').replace(']','').replace("'",""), str(dict_koid[key]).strip('[]').replace(']','').replace('[','').replace("'",""), len(values), sep ='\t', file = fout)
