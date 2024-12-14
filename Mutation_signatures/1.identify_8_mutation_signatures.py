import pandas as pd
import collections 
import argparse
import os
import glob
import itertools
import matplotlib.pyplot as plt
import numpy as np
# from sklearn.decomposition import NMF
from tqdm import tqdm
from sklearn.decomposition import NMF
import random


# -------------------------step1 get the matrix of 96 mutation types-------------------------------------
'''
The format of output of each sample
codons,-,AAA,AAC,AAG,AAT,ACA,ACC,ACG,ACT,AGA,AGC,AGG,AGT,ATA,ATC,ATG,ATT,CAA,CAC,CAG,CAT,CCA,CCC,CCG,CCT,CGA,CGC,CGG,CGT,CTA,CTC,CTG,CTT,GAA,GAC,GAG,GAT,GCA,GCC,GCG,GCT,GGA,GGC,GGG,GGT,GTA,GTC,GTG,GTT,TAA,TAC,TAG,TAT,TCA,TCC,TCG,TCT,TGA,TGC,TGG,TGT,TTA,TTC,TTG,TTT
-,0,1,3,1,2,1,1,3,2,1,0,0,1,2,2,1,1,0,0,0,2,3,0,0,1,0,1,0,0,1,1,1,0,5,3,3,8,7,1,0,0,5,3,2,1,0,0,0,1,0,0,1,2,1,0,0,0,1,1,1,0,0,1,2,2
AAA,7,1,14,908,24,12,0,0,0,91,0,0,0,5,0,0,0,25,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,106,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
AAC,1,22,0,10,796,0,16,0,0,0,67,0,0,0,5,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,98,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0
AAG,2,961,19,0,17,0,0,12,0,0,0,55,0,0,0,4,0,0,0,19,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,31,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0
AAT,3,32,846,13,0,0,1,0,10,0,0,0,94,0,0,0,4,0,0,0,21,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,135,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,11,0,0,0,0,0,0,0,0,0,0,0,0
ACA,2,9,0,0,0,0,95,506,161,1,0,0,0,54,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,113,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,0,0,0,0,0,0,0,0,0,0,0
ACC,1,0,13,0,0,105,0,72,487,0,15,0,0,0,38,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,107,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,11,0,0,0,0,0,0,0,0,0,0
ACG,1,0,0,10,0,469,56,0,91,0,0,3,0,0,0,42,0,0,0,0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,54,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0
ACT,0,0,0,0,11,160,534,94,1,0,0,0,17,0,0,0,44,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,118,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,13,0,0,0,0,0,0,0,0
AGA,2,89,0,0,0,1,0,0,0,0,5,111,5,3,0,0,0,0,0,0,0,0,0,0,0,23,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,19,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2,0,0,0,0,0,0,0
AGC,0,0,94,0,0,0,19,0,0,4,2,4,507,0,1,0,0,0,0,0,0,0,0,0,0,0,4,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,42,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6,0,0,0,0,0,0
...
...


'''
def get_path_of_suffix(main_path,suffix):
    path_list=os.listdir(main_path)
    list=[]
    # for path in path_list:
    #     if os.path.splitext(path)[1]==f'.{suffix}':
    #         list.append(f'{main_path}/{path}')
    list=glob.glob(f'{main_path}/*.{suffix}')
    print(list)
    return list

def remove_N(data):
    df=pd.DataFrame()
    index_list=[]
    for i in data.index:
        name=data.iloc[i,0]
        if 'N' in name:
            index_list.append(i)
        
    if len(index_list) != 0:
        df=data.drop(index=index_list,axis=0)
    else:
        df=data
    df1=pd.DataFrame()
    col_list=[]
    for name in df.columns:
        if 'N' in name:
            print(name)
            col_list.append(name)
    if len(index_list) != 0:
        df1=df.drop(columns=col_list)
    else:
        df1=df 
            
    return df1

def transform_SBS_from_csv(path,mainwritepath):
    pathlist=get_path_of_suffix(path,'csv')
    for path in pathlist:
        name=os.path.basename(path).split('.')[0]
        with open(path,'r') as f:
            lines=f.readlines()
            list_name=[]
            for i in range(len(lines)) :
                line=lines[i]
                line=line.strip()
                
                if line.startswith('codons'):
                    #print(line)
                    list_colums=line.split(',')
                    for n in list_colums:
                        n=n.strip()
                        list_name.append(n)
            
                   
                    #print(list_colums)
                    j=i
                    break
            dict=collections.defaultdict(list)
            for m in range(j+1,len(lines)):
                line=lines[m]
                line=line.strip()
                if line.startswith('# Amino'):
                    break
                if line.strip() == "":
                    break
                else:
                    data=line.split(',')
                    print(len(data))
                    print(len(list_colums))
                    for i in range(len(list_colums)):
                        dict[list_name[i]].append(data[i].strip()) 
                
            dic=pd.DataFrame(dict)
            
            dic=remove_N(dic)
            writrpath=f'{mainwritepath}/{name}.csv'
            dic.to_csv(writrpath,sep=",",index=0)
            #print(dict)
            print(list_name)
            print(list_colums)

transform_SBS_from_csv('/media/desk16/szz077/team/metagene/cancer_snpeff/MEL',r'/media/desk16/szz077/team/metagene/result/96mutation/eff_output')
transform_SBS_from_csv('/media/desk16/szz077/team/metagene/cancer_snpeff/NSCLC',r'/media/desk16/szz077/team/metagene/result/96mutation/eff_output')
transform_SBS_from_csv('/media/desk16/szz077/team/metagene/cancer_snpeff/RCC',r'/media/desk16/szz077/team/metagene/result/96mutation/eff_output')
transform_SBS_from_csv('/media/desk16/szz077/team/metagene/cancer_snpeff/GIC',r'/media/desk16/szz077/team/metagene/result/96mutation/eff_output')


#-----------------------step2 get the matrix of feature of each sample----------------------------
'''
the example of the output of each sample

,C>A,C>G,C>T,T>A,T>C,T>G
ANA,12,2,143,8,132,14
ANC,14,34,132,7,121,18
ANG,11,4,95,14,109,13
ANT,13,31,131,11,127,10
CNA,5,1,30,14,28,3
CNC,5,4,43,2,49,4
CNG,11,2,92,43,83,14
CNT,7,4,54,6,56,5
GNA,28,15,171,7,168,31
GNC,17,20,118,4,120,14
GNG,11,5,90,8,91,10
GNT,12,28,119,7,112,13
TNA,1,0,30,3,24,0
TNC,7,5,26,14,23,5
TNG,11,3,24,1,20,8
TNT,10,2,37,45,37,16
'''
def count_mutation_feature(data):
    dict_AA=collections.defaultdict(list)
    list_all=[f"{i}{j}{n}"  for i,j,n in itertools.product('ACGT',repeat=3)]
    for i in list_all:
        if i not in data.columns:
            data[i]=0
           
    for i in list_all:
        if i not in data.index:
            data.loc[i,:]=0
          

    for i,j in itertools.product('ACGT',repeat=2):
        list_AA=[f'{i}{b}{j}' for b in "ACGT"]

        df=data.loc[list_AA,list_AA]
        list_mutation_mode=['C>A','C>G','C>T','T>A','T>C','T>G']
        dict_AA['C>A'].append(df.iloc[1,0]+df.iloc[2,3])
        '''print(df)
        print(df.iloc[1,0],df.iloc[2,3])
        print(dict_AA)
        break'''
        dict_AA['C>G'].append(df.iloc[1,2] + df.iloc[2,1])
        dict_AA['C>T'].append(df.iloc[1,3] + df.iloc[2,0])
        dict_AA['T>A'].append(df.iloc[3,0] + df.iloc[0,3])
        dict_AA['T>C'].append(df.iloc[3,1] + df.iloc[0,2])
        dict_AA['T>G'].append(df.iloc[3,2] + df.iloc[0,1])

    df_AA=pd.DataFrame(dict_AA)
    index_list=[]
    for i,j in itertools.product('ACGT',repeat=2):
        index_list.append(f'{i}N{j}')
    df_AA.index=index_list
    return df_AA

# processed_data=count_mutation_feature(data)
# print(processed_data.shape)

# get the input martix
path=r'/media/desk16/szz077/team/metagene/result/96mutation/eff_output'
paths=os.listdir(path)
pathlist=[f'{path}/{i}' for i in paths]

for t in tqdm(pathlist):
    try:
        name=os.path.basename(t).split('.')[0]
        print(name)
        data=pd.read_csv(t,sep=',',header=0)
        data.index=data.iloc[:,0]
        processed_data=count_mutation_feature(data)
        processed_data.to_csv(f'/media/desk16/szz077/team/metagene/result/96mutation/feature/{name}.csv',sep=',')
    except Exception as e:
            print(e, type(e))
            if (isinstance(e, pd.errors.EmptyDataError)):
                print("The file is empty")

#-----------------------srep3----------------------------
'''
generated a 96 mutation count matrix V (96*m),where m represents the number of samples
'''
path=r'/media/desk16/szz077/team/metagene/result/96mutation/feature'
paths=os.listdir(path)
pathlist=[f'{path}/{i}' for i in paths]
np_list=[]
col_list=[]
for t in tqdm(pathlist):

    name=os.path.basename(t).split('.')[0]
    col_list.append(name)
    # print(name)
    data=pd.read_csv(t,index_col=0)

    data=data.values
    data=data.reshape(96,-1)
    data=data/(data.sum())
    
    
    np_list.append(data)
a=np.concatenate(np_list,axis=1)
a=np.nan_to_num(a)


#-----------------------step4----------------------------
'''
Using Non-negative Matrix Factorization (NMF), 8 mutation signatures were identified.
'''
# R1=100000000
R_list=[]
model = NMF(n_components=8, init='random', random_state=583)
W = model.fit_transform(a)
H = model.components_
R=model.reconstruction_err_
# print(W)
# print(R)
H=pd.DataFrame(H,columns=col_list)
# print(W)
H.to_csv(r'/media/desk16/szz077/team/metagene/result/96mutation/matrix/H.csv',sep=',',index=None)
np.savetxt('/media/desk16/szz077/team/metagene/result/96mutation/matrix/W.csv', W, delimiter=",")


#-----------------------step5---------------------------------
'''
identify the mutation signature for each patient
'''
data=pd.read_csv(r'/media/desk16/szz077/team/metagene/result/96mutation/matrix/H.csv',sep=',')
meta=pd.read_csv('/media/desk16/szz077/team/metagene/metadate_new.csv',sep=',')
sig_index=[f'S{i}' for i in range(1,9)]
sig_index
data.index=sig_index
data

max_index=data.idxmax()
max_index=max_index.reset_index().rename(columns={'index':'run',0:'sig'})
max_index
df=pd.merge(max_index,meta,on='run',how='left')
df.to_csv('/media/desk16/szz077/team/metagene/result/96mutation/single_sample/max_sig.csv',index=None)