# -*- coding: utf-8 -*-
# @Author: Xue Chao
# @Time: 2025/03/11 10:59
# @Function:
import gzip
import math
import os

import filetype
import numpy as np
import pandas as pd
from scipy.stats import norm

from paper.para import GWAS_raw_meta, GWAS_raw_data_dir, Full_fine_GWAS_dir, hg19_rsid_path, hg38_rsid_path, \
    raw_hg19_rsid_path, raw_hg38_rsid_path, MAGMA_GWAS_dir
from paper.util import log, make_dir


def get_skip_row_count(file_path,comment):
    isGzip = False
    try:
        if str(filetype.guess(file_path).extension) == "gz":
            isGzip = True
    except:
        isGzip = False
    n_rows=0
    if isGzip:
        with gzip.open(file_path, "r") as reader:
            while True:
                line = reader.readline().decode().strip('\n')
                if not line.startswith(comment):
                    break
                n_rows += 1
    else:
        with open(file_path, "r") as reader:
            while True:
                line = reader.readline().strip('\n')
                if not line.startswith(comment):
                    break
                n_rows+=1
    return n_rows

def save_clean_rsid():
    log(f'start')
    rsid_dbs=[{'raw':raw_hg19_rsid_path,'fine':hg19_rsid_path},{'raw':raw_hg38_rsid_path,'fine':hg38_rsid_path}]
    for db in rsid_dbs:
        raw_path=db['raw']
        fine_path=db['fine']
        skip_row=get_skip_row_count(raw_path,comment='##')
        df=pd.read_table(raw_path,skiprows=skip_row,usecols=[0,1,2],dtype={0: 'str', 1: 'int', 2: 'str'})
        print(df)
        df.to_csv(fine_path,sep='\t',index=False,lineterminator='\n')
        log(f'load {df.shape[0]} snps in db.')


def read_rsid_db(ref_geno):
    rsid_path = ''
    if ref_geno == 'hg19':
        rsid_path = hg19_rsid_path
    if ref_geno == 'hg38':
        rsid_path = hg38_rsid_path
    if rsid_path == '':
        raise Exception(f'unsupported ref genome version {ref_geno}')
    df=pd.read_table(rsid_path,dtype={0: 'str', 1: 'int64', 2: 'str'})
    df.rename(columns={'#CHROM':'chr','POS':'bp','ID':'snp'},inplace=True)
    return df


def format_gwas():
    ## GWAS must have chr,bp,p columns
    # cols must have
    make_dir(Full_fine_GWAS_dir)
    must_cols=['snp','chr','bp','a1','a2','z','p','n_eff']
    meta_df=pd.read_excel(GWAS_raw_meta)
    for i in meta_df.index:
        abbr=meta_df.loc[i,'abbr']
        category=meta_df.loc[i,'category']
        fdir=meta_df.loc[i,'dir']
        path=meta_df.loc[i,'file_name']
        sep=meta_df.loc[i,'sep']
        ref_geno=meta_df.loc[i,'ref_genome']
        file_path=f'{GWAS_raw_data_dir}/{category}/{fdir}/{path}'
        skip_row=0
        comment= meta_df.loc[i, 'comment_prefix']
        if not pd.isna(comment) and str(comment).strip()!='':
            skip_row=get_skip_row_count(file_path,str(comment).strip())
        log(f'{abbr}: skip {skip_row} rows')
        df=pd.read_table(file_path,skiprows=skip_row,sep=sep)
        # change normal name
        rename_cols={}
        for c in must_cols:
            cc=meta_df.loc[i,c]
            if not pd.isna(cc) and str(cc).strip() != '':
                rename_cols[str(cc).strip()]=c
        df.rename(columns=rename_cols,inplace=True)
        # add snp columns
        snp_col=meta_df.loc[i,'snp']
        if pd.isna(snp_col) or str(snp_col).strip()=='':
            rsid_df=read_rsid_db(ref_geno)
            log(f'{abbr} [add rsid]: load {rsid_df.shape[0]} SNPs')
            # rsid_mapping = {}
            # for index, row in rsid_df.iterrows():
            #     key = (row['#CHROM'], row['POS'])
            #     rsid_mapping[key] = row['ID']
            # log(f'{abbr}: trans db to dict')
            # def get_rsid(row):
            #     key = (row['chr'], row['bp'])
            #     return rsid_mapping.get(key, f"{row['chr']}:{row['bp']}")
            # df['snp'] = df.apply(get_rsid, axis=1)
            df['chr']=df['chr'].astype('str')
            df['bp']=df['bp'].astype('int64')
            df = pd.merge(df, rsid_df[['chr', 'bp', 'snp']], on=['chr', 'bp'], how='left')
            log(f'{abbr}: merged')
        df['snp'] = df['snp'].fillna(df['chr'].astype(str) + ':' + df['bp'].astype(str))
        log(f'{abbr}: fill SNP nan')
        # add z columns
        z_col=meta_df.loc[i,'z']
        if pd.isna(z_col) or str(z_col).strip()=='':
            beta_col=meta_df.loc[i,'beta']
            if not pd.isna(beta_col) and str(beta_col).strip() != '':
                df.rename(columns={str(beta_col).strip(): 'beta'},inplace=True)
            else:
                or_col = meta_df.loc[i, 'or']
                if pd.isna(or_col) or str(or_col).strip() == '':
                    raise Exception('No beta and no or')
                df.rename(columns={str(or_col).strip(): 'or'},inplace=True)
                df['beta']=np.log(df['or'])
            df['z'] = abs(norm.ppf(1-df['p'] / 2))*np.sign(df['beta'])
            df['se'] = df['beta'] / df['z']

        # add n sample
        n_col = meta_df.loc[i, 'n_eff']
        if pd.isna(n_col) or str(n_col).strip() == '':
            nsam_col = meta_df.loc[i, 'n_sample']
            if pd.isna(nsam_col) or str(nsam_col).strip() == '':
                raise Exception('No sample size')
            df['n_eff']=int(nsam_col)
        all_cols=must_cols.copy()
        for avail_col in ['beta','se']:
            if avail_col in df.columns:
                all_cols.append(avail_col)
        df=df[all_cols]
        df['a1']=df['a1'].str.upper()
        df['a2'] = df['a2'].str.upper()
        print(df)
        df.to_csv(f'{Full_fine_GWAS_dir}/{abbr}.gwas.sum.tsv.gz',index=False,sep='\t',lineterminator='\n',float_format='%.4g',na_rep='NaN')
        log(f'save {abbr}')


def format_gwas_magma():
    wdir=Full_fine_GWAS_dir
    odir = MAGMA_GWAS_dir
    make_dir(odir)
    for f in os.listdir(wdir):
        if not f.endswith('.tsv.gz'):
            continue
        log(f'start {f.split(".")[0]}')
        df=pd.read_table(f'{wdir}/{f}',dtype=str)
        df=df[['snp','chr','bp','p','n_eff']]
        df.columns=['SNP','CHR','BP','P','NOBS']
        out_path=f'{odir}/{f.split(".gz")[0]}'
        df.to_csv(out_path,sep='\t',index=False,lineterminator='\n')
        log(f'save to {out_path}')


def main():
    # save_clean_rsid()
    # read_rsid_db('hg19')
    # format_gwas()
    format_gwas_magma()

if __name__ == '__main__':
    main()