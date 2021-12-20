#!/usr/bin/env python
# coding: utf-8

# In[1]:
import importlib.util, sys

def findDeps(name):
    if (spec := importlib.util.find_spec(name)) is not None: return
    else:
        print(f"""
            can't find the {name!r} module
            you can manually install the dependencies using:

                pip install {name}

                or 

                conda install {name}
        """)
        sys.exit()


for i in ['pandas', 'numpy']: findDeps(i)
import pandas as pd
import numpy as np

vcfIn, gtfIn, saveFile = "", "", 'report.txt'
def validateExtension(fil):
    if fil.find('.gz') > 0:
        print("""
        \033[1;31mERROR: This program does not accept '.gz' files%\033[0m
        Please input uncompressed files

        use -h or --help for more information
        """)
        sys.exit()
    else: return


if len(sys.argv)==4: 
    vcfIn = sys.argv[1]
    gtfIn = sys.argv[2]
    saveFile = sys.argv[2]
    validateExtension(vcfIn)
    validateExtension(gtfIn)

elif len(sys.argv)==3: 
    vcfIn = sys.argv[1]
    gtfIn = sys.argv[2]
    validateExtension(vcfIn)
    validateExtension(gtfIn)
elif len(sys.argv)==2:
    if sys.argv[1] == "-h" or sys.argv[1] == "--help":
        print("""
        --------------------------------------------------------------------
        vcfAnotator
        --------------------------------------------------------------------\n
        This script is built with the intention of annotating VCF files
        inorder to identify the genes that the identified mutations affect.

        command line arguments accepted:   
        1. vcf_file_name
        2. gtf_file_name
        3. output_file_name (Optional, default out put file is 'report.txt')

        or 

        use -h or --help for more information
        --------------------------------------------------------------------
        """)
        
        sys.exit()
else: 
    print("""
    \033[1;31mERROR: input files not provided as arguments%\033[0m
    please enter the file names in the given order:
        1. vcf_file_name
        2. gtf_file_name
        3. output_file_name (Optional, default out put file is 'report.txt')

        use -h or --help for more information
    """)
    sys.exit()

# In[2]:


def getStart(loc):
    with open(loc) as f:
        lines, fr = f.readlines(), 0        
        for i in lines:
            fr +=1
            if i[0]=="#": continue
            else: break
        return fr-1
        
def procDesc(inDf, srch, filt=True, validator="protein_coding"):
    ret = []
    desc = list(inDf['desc'])
    for i in desc: 
        proc = i[:-1]
        obs = [j.split() for j in proc.split(';')]
        obs = np.array(obs,dtype=object).T.tolist()
        keylist = dict(zip(obs[0], obs[1]))
        
        #if set to filter mode the arrat is populated with booleans
        if filt:
            if 'gene_biotype' in obs[0]: 
                if keylist['gene_biotype'].find(validator)<0: ret.append(1)
                else: ret.append(0)
            else: ret.append(0)
                
        # else the array is populated with the key values for the search term
        else: 
            if srch in keylist:
                ret.append(keylist[srch])
            else: ret.append('N/A')
    return ret


# In[3]:


gtfIn = './data/Schizosaccharomyces_pombe.ASM294v2.51.gtf'

nams = ['chrom', 'source', 'type', 'strt', 'end', '1', 'dir', '3', 'desc']
# importing gff as gff
gff = pd.read_csv(gtfIn, sep="\t", skiprows=getStart(gtfIn), names=nams)

#filtering out non gene info
gff = gff[gff['type']=="gene"]
gff['protcod'] = procDesc(gff,'protein_coding')
gff = gff[gff['protcod']==1].reset_index(drop=True) # filtering out non protein coding genes

gff['gene_id'] = procDesc(gff,'gene_id', False)
gff['gene_name'] = procDesc(gff,'gene_name', False)


# In[5]:


# process vcf file
def getOneLine(fil, pos): 
    ret = "" 
    with open(fil) as f: 
        for i in range(pos+1):
            if i != pos: ret = f.readline()
            else: break
    return ret

def getVcf(fil):
    nams = [i.lower() for i in getOneLine(fil, getStart(fil))[1:-1].split('\t')]
    return pd.read_csv(fil, sep="\t", skiprows=getStart(fil), names=nams)

vcf = getVcf(vcfIn)
vcf = vcf[['chrom', 'pos']]

vcfArr = vcf.to_numpy()


# In[24]:


ids = []
nams = []
for i in vcfArr: 
    cr, po = i
    res = gff[gff['chrom']==cr]
    res = gff[gff['strt']<=po]
    res = gff[gff['end']>=po]
    ids.append(list(res['gene_id']))
    nams.append(list(res['gene_name'])) 

#remove 'N/A' alues
nams=list(map(lambda x: list(filter(lambda a: a!='N/A', x)), nams))

vcf['gene_id'] = list(map(lambda x: ';'.join(list(map(lambda a: a.replace('"', '') , x))), ids))
vcf['gene_name'] =list(map(lambda x: ';'.join(list(map(lambda a: a.replace('"', '') , x))), nams))
vcf = vcf.reset_index(drop=True)


# In[25]:


vcf.to_csv(saveFile, sep="\t")


# %%
