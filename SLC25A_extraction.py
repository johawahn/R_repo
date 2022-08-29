#/Users/johannawahn/Desktop/iGEM/Modeling
# -*- coding: utf-8 -*-
"""
Created on %(29/6/2020)s

@author: %(johananWahn)s
"""

import re
import os
import ntpath
import pandas as pd
import csv



#############################################################################################################
SLC25A39 = "ENSG00000013306"
SLC25A40 = "ENSG00000075303"


def extraction(dir, prot):      
    #Change directory
    os.chdir(dir)
    
    
    
    #Creating the data frame
    dim = list(range(0,1000)) #Number max of samples
    dict = {"Sample": dim, "Protein": dim, "FPKM": dim}
    df = pd.DataFrame(dict)
    index = 0
    nb_file = 0

    
    for file in os.listdir(dir):
        
        if ntpath.basename(file) == ".DS_Store":
            continue
        
        else:
            with open(file) as fl:
                for line in fl:
                    if re.search(prot, line):
                        expr = line.split()
                        
                        df['Sample'][index] = ntpath.basename(file)
                        df['Protein'][index] = expr[0]
                        df['FPKM'][index] = float(expr[1])
                        
                        index = index + 1
                        nb_file = nb_file + 1
    
    df = df.loc[range(0,nb_file),:] 
    return(df)


All_FPKM_UQ_COAD_SLC25A39 = extraction("/Users/johannawahn/Desktop/Zamboni_Lab/Transcriptome_CRC/TXT_FPKM_UQ_COAD", SLC25A39)                                                 
KRAS_G12C_SLC25A39 = extraction("/Users/johannawahn/Desktop/Zamboni_Lab/Transcriptome_CRC/TXT_KRAS_G12C", SLC25A39)
KRAS_G12V_SLC25A39 = extraction("/Users/johannawahn/Desktop/Zamboni_Lab/Transcriptome_CRC/TXT_KRAS_G12V", SLC25A39)
KRAS_G13D_SLC25A39 = extraction("/Users/johannawahn/Desktop/Zamboni_Lab/Transcriptome_CRC/TXT_KRAS_G13D", SLC25A39)
NORMAL_SLC25A39 = extraction("/Users/johannawahn/Desktop/Zamboni_Lab/Transcriptome_CRC/TXT_NORMAL", SLC25A39)

All_FPKM_UQ_COAD_SLC25A40 = extraction("/Users/johannawahn/Desktop/Zamboni_Lab/Transcriptome_CRC/TXT_FPKM_UQ_COAD", SLC25A40)
KRAS_G12C_SLC25A40 = extraction("/Users/johannawahn/Desktop/Zamboni_Lab/Transcriptome_CRC/TXT_KRAS_G12C", SLC25A40)
KRAS_G12V_SLC25A40 = extraction("/Users/johannawahn/Desktop/Zamboni_Lab/Transcriptome_CRC/TXT_KRAS_G12V", SLC25A40)
KRAS_G13D_SLC25A40 = extraction("/Users/johannawahn/Desktop/Zamboni_Lab/Transcriptome_CRC/TXT_KRAS_G13D", SLC25A40)
NORMAL_SLC25A40 = extraction("/Users/johannawahn/Desktop/Zamboni_Lab/Transcriptome_CRC/TXT_NORMAL", SLC25A40)

####################################### SAVE AS CSV FILE ####################################################

All_FPKM_UQ_COAD_SLC25A39.to_csv('All_FPKM_UQ_COAD_SLC25A39.csv')
KRAS_G12C_SLC25A39.to_csv('KRAS_G12C_SLC25A39.csv')
KRAS_G12V_SLC25A39.to_csv('KRAS_G12V_SLC25A39.csv')
KRAS_G13D_SLC25A39.to_csv('KRAS_G13D_SLC25A39.csv')
NORMAL_SLC25A39.to_csv('NORMAL_SLC25A39.csv')

All_FPKM_UQ_COAD_SLC25A40.to_csv('All_FPKM_UQ_COAD_SLC25A40.csv')
KRAS_G12C_SLC25A40.to_csv('KRAS_G12C_SLC25A40.csv')
KRAS_G12V_SLC25A40.to_csv('KRAS_G12V_SLC25A40.csv')
KRAS_G13D_SLC25A40.to_csv('KRAS_G13D_SLC25A40.csv')
NORMAL_SLC25A40.to_csv('NORMAL_SLC25A40.csv')


############################################### GRAPH ####################################################

import seaborn as sns
import matplotlib.pyplot as plt
sns.set(style="whitegrid")

f, ax = plt.subplots(figsize=(8, 8))


####  SLC25A39
dict = {"ALL": dim, "Protein": dim, "FPKM": dim}
df_SLC25A39 =  = pd.DataFrame(dict)


# Show each distribution with both violins and points
sns.violinplot(x="Cell type",y="FPLM",data=chick, inner="box", palette="Set3", cut=2, linewidth=3)

sns.despine(left=True)

f.suptitle('FPKM by cell type', fontsize=18, fontweight='bold')
ax.set_xlabel("Cell type",size = 16,alpha=0.7)
ax.set_ylabel("FPKM",size = 16,alpha=0.7)





