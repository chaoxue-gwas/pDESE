# -*- coding: utf-8 -*-
# @Author: Xue Chao
# @Time: 2025/03/12 10:16
# @Function:
import os

import pandas as pd

from paper.para import output_DIR, magma_gene_annot, magma_ref_geno, MAGMA_GWAS_dir, pDESE_combine_data_noNA_DIR, \
    pDESE_combine_data_magma_DIR, pDESE_combine_data_DIR, pDESE_combine_rez_DIR
from paper.util import make_dir, LOCAL_DIR, batch_shell_plus, run_command, log, extract_col

magma_bin=f'{LOCAL_DIR}/lib/magma/magma'


class MAGMA:
    def __init__(self,nt):
        self.nt=nt
        self.gwas_suffix='.sum.tsv'
        self.out_dir=f'{output_DIR}/magma'
        make_dir(self.out_dir)

    def magma_gene(self):
        cmds=[]
        for gwas in os.listdir(MAGMA_GWAS_dir):
            if not gwas.endswith(self.gwas_suffix):
                continue
            pheno=gwas.split('.')[0]
            magma_gwas_path=f'{MAGMA_GWAS_dir}/{gwas}'
            result_prefix=f'{self.out_dir}/{pheno}'
            sys_log = result_prefix + '.sys_log'
            make_dir(os.path.dirname(result_prefix))
            cmd1 = f'{magma_bin} --annotate --snp-loc {magma_gwas_path} --gene-loc {magma_gene_annot}' \
                   f' --out {result_prefix}'
            cmd2 = f'{magma_bin} --bfile {magma_ref_geno} --gene-annot {result_prefix}.genes.annot ' \
                   f'--pval {magma_gwas_path} ncol=NOBS --out {result_prefix}'
            cmd = f'{cmd1} && {cmd2}'
            # run_command(cmd, sys_log)
            cmds.append(cmd)
        batch_shell_plus(cmds,self.nt)


    def magma_tissue(self):
        edir=pDESE_combine_data_magma_DIR
        expr_files=[]
        for f in os.listdir(edir):
            expr_files.append({'name':'.'.join(f.split('.')[:2]),'path':f'{edir}/{f}'})
        cmds=[]
        for gwas in os.listdir(MAGMA_GWAS_dir):
            if not gwas.endswith(self.gwas_suffix):
                continue
            pheno=gwas.split('.')[0]
            result_prefix=f'{self.out_dir}/{pheno}'
            for expr in expr_files:
                cmd = f'{magma_bin} --gene-results {result_prefix}.genes.raw ' \
                       f'--gene-covar {expr["path"]} --out {result_prefix}.{expr["name"]}'
                cmds.append(cmd)
        batch_shell_plus(cmds,self.nt)

    def magma_gene_change(self):
        gdf=pd.read_table(magma_gene_annot,header=None,dtype=str)
        id_gene={}
        for i in gdf.index:
            id_gene[gdf.loc[i,0]]=gdf.loc[i,5]
        for gwas in os.listdir(MAGMA_GWAS_dir):
            if not gwas.endswith(self.gwas_suffix):
                continue
            pheno=gwas.split('.')[0]
            result_prefix=f'{self.out_dir}/{pheno}'
            gr=f'{result_prefix}.genes.out'
            gn=f'{result_prefix}.genes.out.symbol'
            df=pd.read_table(gr,dtype=str,sep='\s+')
            print(df)
            print(df.columns)
            df['GENE']=df['GENE'].map(lambda x:id_gene[x])
            df.to_csv(gn,sep='\t',index=False,lineterminator='\n')
            log(f'save {pheno}')

    def magma_tissue_expr(self):
        odir = pDESE_combine_data_magma_DIR
        make_dir(odir)
        gdf=pd.read_table(magma_gene_annot,header=None,dtype=str)
        id_gene={}
        for i in gdf.index:
            id_gene[gdf.loc[i,5]]=gdf.loc[i,0]
        wdir = pDESE_combine_rez_DIR
        for f in os.listdir(wdir):
            if not f.endswith('.txt.gz'):
                continue
            log(f'{f} start')
            raw_df = pd.read_table(f'{wdir}/{f}', index_col=0)
            df = extract_col(raw_df, '.z')
            df.index=df.index.map(lambda x:id_gene.get(x,x))
            unique_index = ~df.index.duplicated()
            df = df[unique_index]
            # df.drop_duplicates(subset=['Gene'],inplace=True)
            df.index.name='Gene'
            df.to_csv(f'{odir}/{f.split(".gz")[0]}', na_rep='0', sep='\t', lineterminator='\n')


def main():
    magma = MAGMA(10)
    # magma.magma_gene()
    # magma.magma_gene_change()
    # magma.magma_tissue_expr()
    magma.magma_tissue()


if __name__ == '__main__':
    main()
