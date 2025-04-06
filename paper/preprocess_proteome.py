# -*- coding: utf-8 -*-
# @Author  : Xue Chao
# @Time    : 2024/03/13 19:28
# @Function: Preprocess Proteome from GTEx, HPA, et al.
import os.path
import re

import pandas as pd

from para import RAW_HPA_PATH, MSB_2019_PRO_PATH, MSB_2019_TRA_PATH, RAW_MSB_2019_PRO_PATH, RAW_MSB_2019_TRA_PATH, \
    RAW_CELL_2020_PRO_PATH, RAW_CELL_2020_TRA_PATH, CELL_2020_PRO_PATH, CELL_2020_TRA_PATH, pDESE_combine_data_DIR, \
    pDESE_combine_data_noNA_DIR, pDESE_combine_data_magma_DIR
from paper.util import make_dir, log, GeneID


def hpa_proteome():
    '''
    raw data: https://www.proteinatlas.org/download/normal_tissue.tsv.zip
    processed data format: matrix and code expression level to number.
    :return:
    '''
    df=pd.read_table(RAW_HPA_PATH)
    df=df.loc[df['Reliability']!='Uncertain',:]

    # for c in ['Tissue', 'Cell type', 'Level', 'Reliability']:
    #     print(df[c].unique())

    pass  

def msb_2019_proteome_and_transcriptome():
    '''
    raw data in is paper supp tables. https://www.embopress.org/doi/full/10.15252/msb.20188503#msb188503-sup-0003
    :return:
    '''
    pro_df=pd.read_table(RAW_MSB_2019_PRO_PATH) # format(tab): Gene ID,Gene name,tissue1-n...
    tra_df=pd.read_table(RAW_MSB_2019_TRA_PATH) # format(tab): Gene ID,Gene name,tissue1-n...

    batch=[{'data':pro_df,'out':MSB_2019_PRO_PATH},{'data':tra_df,'out':MSB_2019_TRA_PATH}]
    for b in batch:
        df,out_path=b['data'],b['out']
        make_dir(os.path.dirname(out_path))
        df.index=df['Gene name']
        df.index.name='Gene'
        df=df[[c for c in df.columns if c not in ['Gene ID','Gene name']]]
        df=df[df.columns[:-4]]
        df.to_csv(out_path,sep='\t',na_rep='nan')

def cell_2020_proteome_and_transcriptome():
    '''
    raw data in is paper supp tables. https://www.cell.com/cell/fulltext/S0092-8674(20)31078-3#secsectitle0095
    :return:
    '''
    pro_df=pd.read_excel(RAW_CELL_2020_PRO_PATH) # format(tab): gene.id,tissue1-n...
    tra_df=pd.read_excel(RAW_CELL_2020_TRA_PATH) # format(tab): gene.id,tissue1-n...
    batch=[{'data':pro_df,'out':CELL_2020_PRO_PATH},{'data':tra_df,'out':CELL_2020_TRA_PATH}]
    for b in batch:
        df,out_path=b['data'],b['out']
        make_dir(os.path.dirname(out_path))
        df.index=df['gene.id']
        df.index.name='Gene'
        df=df[[c for c in df.columns if c not in ['gene.id']]]
        df.to_csv(out_path,sep='\t',na_rep='nan')


def norm_pDESE_data():
    wdir=pDESE_combine_data_DIR
    make_dir(wdir)
    wpaths=[{'db_name':'MSB_2019','pro':MSB_2019_PRO_PATH,'tra':MSB_2019_TRA_PATH},
            {'db_name':'CELL_2020','pro':CELL_2020_PRO_PATH,'tra':CELL_2020_TRA_PATH}]
    for it in wpaths:
        db=it['db_name']
        pro=it['pro']
        tra=it['tra']
        pro_df=pd.read_table(pro,index_col=0)
        tra_df=pd.read_table(tra,index_col=0)
        pro_df=pro_df[~pro_df.index.duplicated()]
        tra_df=tra_df[~tra_df.index.duplicated()]
        common_gene=sorted(set(pro_df.index).intersection(tra_df.index))
        gid=GeneID()
        gene_symbol=gid.get_gene_symbol(common_gene)
        eff_gs=[]
        eff_cg=[]
        for i in range(len(common_gene)):
            if gene_symbol[i]!='NA':
                eff_cg.append(common_gene[i])
                eff_gs.append(gene_symbol[i])
        log(f'pro gene: {pro_df.shape[0]}; tra gene: {tra_df.shape[0]}; common gene: {len(common_gene)}; symbol: {len(eff_cg)}')
        pro_df=pro_df.loc[eff_cg,:]
        tra_df=tra_df.loc[eff_cg,:]
        pro_df.index=eff_gs
        tra_df.index=eff_gs
        pro_df.index.name='Gene'
        tra_df.index.name='Gene'
        pro_df.columns=pro_df.columns.map(lambda x:re.sub(r"[\\/\s+]", "", x))
        tra_df.columns = tra_df.columns.map(lambda x: re.sub(r"[\\/\s+]", "", x))
        pro_df.to_csv(f'{wdir}/{db}.proteome.txt.gz',sep='\t',na_rep='NaN')
        tra_df.to_csv(f'{wdir}/{db}.transcriptome.txt.gz',sep='\t',na_rep='NaN')
        log(f'save {db}')


def norm_pDESE_noNA():
    indir = pDESE_combine_data_DIR
    odir = pDESE_combine_data_noNA_DIR
    make_dir(odir)
    for f in os.listdir(indir):
        df=pd.read_table(f'{indir}/{f}',dtype=str)
        df.to_csv(f'{odir}/{f}',na_rep='0',sep='\t',lineterminator='\n',index=False)

def norm_MAGMA():
    indir = pDESE_combine_data_DIR
    odir = pDESE_combine_data_magma_DIR
    make_dir(odir)
    for f in os.listdir(indir):
        df=pd.read_table(f'{indir}/{f}',dtype=str)
        df.to_csv(f'{odir}/{f.split(".gz")[0]}',na_rep='0',sep='\t',lineterminator='\n',index=False)

def main():
    # hpa_proteome()
    # msb_2019_proteome_and_transcriptome()
    # cell_2020_proteome_and_transcriptome()
    # norm_pDESE_data()
    # norm_pDESE_noNA()
    norm_MAGMA()

if __name__ == '__main__':
    main()
