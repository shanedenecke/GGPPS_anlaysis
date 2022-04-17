#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 24 08:53:34 2020

@author: shanedenecke
"""

import os
import pandas as pd
from Bio import SeqIO 
import re
import paramiko
from scp import SCPClient   


os.chdir('/home/sdenecke/Dropbox/Shane/target_discovery/Lep_Target_Characterization/GGPPS/OrthoDB_analysis')


odb_fa=SeqIO.to_dict(SeqIO.parse('/home/sdenecke/Dropbox/Shane/Applications/Custom_Applications/OrthoDB_source/odb_arthropoda_augment.faa',format='fasta'))
xrefs=pd.read_csv('/home/sdenecke/Dropbox/Shane/Applications/Custom_Applications/OrthoDB_source/odb10v1_gene_xrefs_extended_arth_sub.tab',sep='\t',header=None)
ogs=pd.read_csv('/home/sdenecke/Dropbox/Shane/Applications/Custom_Applications/OrthoDB_source/odb10v0_OG2genes.33208.tab',sep='\t',header=None)
taxidKey=pd.read_csv('./inputs/taxid_sp_convert.tsv',sep='\t',header=None)
supp_fa=list(SeqIO.parse('./inputs/Supp_GGPPS_sequences.faa',format='fasta'))
species_key=pd.read_csv('./inputs/species_key.csv')

with open('./inputs/Lepidoptera_species.txt') as f:
    leps = [line.rstrip() for line in f]

with open('./inputs/Pest_species.txt') as f:
    pests = [line.rstrip() for line in f]

with open('./inputs/Pollinator_species.txt') as f:
    pols = [line.rstrip() for line in f]
  
all_species=list(set(pols+pests+leps))
all_species.remove('Chilo_suppressalis')
allspecies_search='|'.join(all_species)


### import drosophila names
dm_ggpps='FBgn0019662'

### subset out orthogroups
ggpps_og=ogs[ogs[1]==xrefs[xrefs[1]==dm_ggpps].iloc[0,0]][0].item()


#### get gene list for each orthogroup
ggpps_genes=list(ogs[ogs[0]==ggpps_og].iloc[:,1])

###CNV
d=dict()
for i in ggpps_genes:
    tax=re.findall('(^.+):.+',i)[0]
    d[tax]=re.findall('^.+:(.+$)',i)[0]

gene_table=pd.DataFrame.from_dict(d,orient='index')
gene_table=gene_table.reset_index()
gene_table.columns=['TaxID','Gene']
gene_table.to_csv('./outputs/GGPPS_raw_genes.csv')


#### Clean up gene table
taxidKey.columns=['TaxID','Abbreviation','Species']
sum_table=pd.merge(gene_table,taxidKey,on='TaxID')
fin_table=sum_table.groupby(['TaxID','Abbreviation','Species'])['Gene'].count()
fin_table.to_csv('./outputs/GGPPS_orthology_table.csv')
merge_table=pd.merge(fin_table,species_key,on='Abbreviation')
merge_table.to_csv('./outputs/OatP74D_Table.csv',index=None)


#### Output table of GGPPS presence in all species
frequency={}
for item in merge_table['Taxonomic Classification']:
   # checking the element in dictionary
   if item in frequency:
      # incrementing the counr
      frequency[item] += 1
   else:
      # initializing the count
      frequency[item] = 1

a=[x for x in frequency.items() if x[1]>10]
b=[x[0] for x in a]

merge_table_sub=merge_table[merge_table['Taxonomic Classification'].isin(b)]



### subset fasta
ggpps_fa=[odb_fa[x] for x in ggpps_genes if x in odb_fa.keys()]

### rename fasta with taxon name
def fasta_rename(seqlist):
    for i in seqlist:
        xref_sub=xrefs[xrefs[0]==i.id]
        if xref_sub.shape[0]>0:
            odb_code=i.id.split(':')[1]
            taxid=i.id.split(':')[0]
            
            species=taxidKey[taxidKey['TaxID']==taxid].iloc[0,2]
            i.id=re.sub(taxid,species,i.id)
            
            if 'ENSEMBL' in set(xref_sub[2]):
                i.id=re.sub(odb_code,xref_sub[xref_sub[2]=='ENSEMBL'].iloc[0,1],i.id) 
            elif 'NCBIgid' in set(xref_sub[2]):
                i.id=re.sub(odb_code,'LOC'+xref_sub[xref_sub[2]=='NCBIgid'].iloc[0,1],i.id) 
            elif 'UniProt' in set(xref_sub[2]):
                i.id=re.sub(odb_code,xref_sub[xref_sub[2]=='UniProt'].iloc[0,1],i.id)     
            elif 'NCBIproteinGI' in set(xref_sub[2]):
                i.id=re.sub(odb_code,xref_sub[xref_sub[2]=='NCBIproteinGI'].iloc[0,1],i.id)     
            elif 'InterPro' in set(xref_sub[2]):
                i.id=re.sub(odb_code,xref_sub[xref_sub[2]=='InterPro'].iloc[0,1],i.id)       
            elif 'Trinity' in set(xref_sub[2]):
                i.id=re.sub(odb_code,xref_sub[xref_sub[2]=='Trinity'].iloc[0,1],i.id)
            else:   
                i.id=i.id
        
    return seqlist
   


     
ggpps_fa_rename=fasta_rename(ggpps_fa)
#npc1b_fa_rename=fasta_rename(npc1b_fa)

ggpps_fa_rename2=ggpps_fa_rename+supp_fa

### Need to add species subset 
ggpps_fa_rename3=[x for x in ggpps_fa_rename2 if re.search(allspecies_search,x.id)]




def fasta_write(l,name):
    with open(name, 'w') as fp:
               for i in l:
                   fp.write('>')
                   fp.write(i.id)
                   fp.write('\n')
                   fp.write(str(i.seq))
                   fp.write('\n')

fasta_write(ggpps_fa_rename2,'./outputs/GGPPS_named_full.faa')
fasta_write(ggpps_fa_rename3,'./outputs/GGPPS_named.faa')
os.system('/home/shanedenecke/Applications/Custom_Applications/fasta_clean.sh -proteome ./outputs/GGPPS_named.faa')
os.system('unset MAFFT_BINARIES; mafft ./outputs/GGPPS_named.faa > ./outputs/GGPPS_named.faa.aln')
os.system('sed -i "s/:/__/g" ./outputs/GGPPS_named.faa.aln')
os.system('~/Applications/trimal/source/trimal -automated1 -fasta -in ./outputs/GGPPS_named.faa.aln -out ./outputs/GGPPS_named.faa.aln.trimm')



# declare credentials   
host = 'chrysalida.imbb.forth.gr'   
username = 'sdenecke'   
password = 'NOT THAT DUMB'
   
# connect to server   
con = paramiko.SSHClient()   
con.load_system_host_keys()   
con.connect(host, username=username, password=password,port=2222) 

with SCPClient(con.get_transport()) as scp:
    scp.put('./outputs/GGPPS_named.faa', '/mnt/disk/shane/temp/')

stdin, stdout, stderr = con.exec_command(' ~/Applications/Custom_Applications/align_and_tree.sh -proteins /mnt/disk/shane/temp/GGPPS_named.faa -threads 10')

for line in stdout.readlines():   
        print(line.strip())
        
with SCPClient(con.get_transport()) as scp:
    scp.get('/home/sdenecke/GGPPS_named/GGPPS_named.nwk.raxml.support','./outputs/GGPPS_tree.nwk')