import pandas as pd
import collections
import os
from tqdm import tqdm
import glob

dict_map = {
    "AT": "SNP:A>T|T>A",
    "TA": "SNP:A>T|T>A",
    "AG": "SNP:A>G|T>C",
    "TC": "SNP:A>G|T>C",
    "AC": "SNP:A>C|T>G",
    "TG": "SNP:A>C|T>G",
    "GA": "SNP:G>A|C>T",
    "CT": "SNP:G>A|C>T",
    "GT": "SNP:G>T|C>A",
    "CA": "SNP:G>T|C>A",
    "GC": "SNP:G>C|C>G",
    "CG": "SNP:G>C|C>G"
}

def get_path_of_suffix(main_path,suffix):
    list=[]
    list=glob.glob(f'{main_path}/*.{suffix}')
    return list

def filter_infr_from_vcf(data):
    type_list=[]
    for i in range(len(data["CHROM"])):
        df=data.iloc[i,:]
        type=f'{df[3]}{df[4]}'
   
        if type not in list(dict_map.keys()):
            type=f'{df[3]}{df[4]}'
        else:
            type=dict_map[f'{df[3]}{df[4]}']
        
        type_list.append(type)
        
    data['alt_type']=type_list
    return data

def count_SNPgroup_number(data):
    SNP_vars=["SNP:A>T|T>A","SNP:A>G|T>C","SNP:A>C|T>G", "SNP:G>A|C>T","SNP:G>T|C>A","SNP:G>C|C>G"]
    # data_R=f"{filemainpathway}/processed_data/common_SNP.csv"
    # data_R=pd.read_csv(data_R,sep=',')
    # R_SNPs=data['alt_type'].value_counts()
    df_R=collections.defaultdict(list)
    #count the munber (Sreies)
    R_number=data["alt_type"].value_counts()

    for i in SNP_vars:
        if i not in R_number.keys():
            df_R[i].append(0)
        else:
            df_R[i].append(R_number[i])
    df_R=pd.DataFrame(df_R)
    return df_R

def count_6muntaitons(path):
    try:
        print(path)
        data=pd.read_csv(path,sep=',')
        data_filter=data[data["INFO"].str.contains('missense_variant')]
        data=filter_infr_from_vcf(data_filter)
        count_data=count_SNPgroup_number(data)
        name=os.path.basename(path)
        name=name.split('.')[0]
        count_data['run']=name

        return count_data
    except Exception as e:
            print(e, type(e))
            if (isinstance(e, pd.errors.EmptyDataError)):
                print("the file is Empty")

    # data_filter=data[data["INFO"].str.contains('missense_variant')]
    # data=filter_infr_from_vcf(data_filter)
    # count_data=count_SNPgroup_number(data)

    # name=os.path.basename(path)
    # name=name.split('.')[0]
    # count_data['run']=name
    # return count_data

path=r'/media/desk16/szz077/team/metagene/result/6mutation/vcf_csv_list'
paths=os.listdir(path)
pathlist=[f'{path}/{i}' for i in paths]

for cancer in pathlist:
    cancer_name=cancer.split('/')[-1]
    path_list=get_path_of_suffix(main_path=f'{cancer}',suffix='csv')

    write_path_file = f'/media/desk16/szz077/team/metagene/result/6mutation/count_result/{cancer_name}_mutation.csv'
    
    df_count=pd.DataFrame()

    # 每个癌症下面的目录
    for path in tqdm(path_list):
        df=count_6muntaitons(path)
        df_count=pd.concat([df_count,df],ignore_index=None)
        print(df_count)
    df_count.to_csv(f'{write_path_file}',index=None,sep=',')
