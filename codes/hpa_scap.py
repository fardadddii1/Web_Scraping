# -*- coding: utf-8 -*-
"""
Created on Thu May 19 11:52:17 2022

@author: basanto
"""

import os
import time
import glob
import numpy as np
import pandas as pd
import requests
from bs4 import BeautifulSoup

#set cwd
wd = ('C:\\Users\\basanto\\Downloads\\histomics_transcriptomics_project\\transcriptomics\\')
os.chdir(wd)

#create data path and dir
data_path = wd+'edgeR_analysis\\0.5_exp\\QLF_Tests\\'
data_dir = data_path+'*.csv'
data_dir = glob.glob(data_dir)

#create general output path
out_path = data_path+'decoded'
if not os.path.isdir(out_path):
    os.mkdir(out_path)

#iterate through data for each study
    
for f in range(len(data_dir)):
    
    file = pd.read_csv(data_dir[f])

    #extract input or scrapping
    gene_tags = file['genes']
    gene_names = file['names']
    num_genes = len(gene_tags)
    
    #extract other expression data
    pval = file['PValue']
    
    #create palceholders
    gene_list = []
    study_data = []
    
    #iterate through genes
    for g in range(len(gene_tags)):
        
        #build search path
        tag = gene_tags[g]
        tag = tag.split('.')
        tag = tag[0]
        name = gene_names[g]
        gene_list.append(name)
        slink = 'https://www.proteinatlas.org/'
        ssearch = slink+tag+'-'+name
        
        #initiate search
        r = requests.get(ssearch)
        
        #create an object to parse the HTML format
        soup = BeautifulSoup(r.content, 'html.parser')
        
        #create placeholder
        gene_data = []
        
        #retrieve all gene details
        entries = soup.find_all('table', class_='border dark round')
        #for entry in entries:
        for entry in range(len(entries)):
            if entry==1:
                continue
            rows = entries[entry].find_all('tr')
            run_tot = 0
            for t in range(1,len(rows)):
                test = (rows[t].find_all('td')[0].text).split('\n')[0]
                gene_data.append(test)       
        study_data.append(gene_data)
        time.sleep(3)
        print('Completed gene no. '+str(g+1)+'/'+str(num_genes))
        
        
    #make gene data uniform
    length = max(map(len,study_data))
    final_data = np.array([xi+[None]*(length-len(xi)) for xi in study_data])  
    
    cols_18 = ['protein','gene_name','tissueSpec','tissueClust','scSpec','scClust',
               'immunSpec','brainSpec','cxProg','Locn','Subcell','protFx','mlclFx',
               'Ds','Detailed_Summary','Etc','Etc2','Etc3']
    cols_17 = cols_18[0:-1] 
    if final_data.shape[1] == 17:
        final_data = pd.DataFrame(final_data,index=gene_list,columns=cols_17)
    else:
        final_data = pd.DataFrame(final_data,index=gene_list,columns=cols_18)
    
    file_out = (data_dir[f].split('\\')[-1].split('.')[0])+'_decoded.xlsx'
    final_data.to_excel(out_path+'\\'+file_out)
    
    print('\n')
    print('\n Completed Study ID: '+(data_dir[1].split('\\')[-1].split('.')[0]))
    print('\n')
       
print('\n')
print('\n All Done !')
print('\n')

