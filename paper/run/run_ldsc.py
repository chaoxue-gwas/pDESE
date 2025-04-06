# -*- coding: utf-8 -*-
# @Author: Xue Chao
# @Time: 2025/02/26 17:51
# @Function:
import os
import re

import pandas
import pandas as pd

from paper.para import pDESE_combine_rez_DIR, output_DIR, gene_annot_ldsc, eur_ref_plink, hg19_gtf, \
    raw_snp_list_ldsc, fine_hm3_snplist, GWAS_raw_meta, GWAS_raw_data_dir, LDSC_GWAS_dir, LDSC_1kg_baseline, \
    LDSC_weights_hm3, Full_fine_GWAS_dir
from paper.util import LOCAL_DIR, extract_col, make_dir, log, run_command, batch_shell_plus

# static parameter
LDSC_SEG_top_gene=1000
LDSC_bin_dir=f'{LOCAL_DIR}/lib/ldsc'



class LDSC:
    def __init__(self,nt:int):
        self.nt=nt
        self.out_dir=f'{output_DIR}/ldsc'
        self.top_gene_dir=f'{self.out_dir}/tissue_top_gene'
        self.tissue_ldsc_dir=f'{self.out_dir}/tissue_ldsc'
        self.ldsc_cts_dir=f'{self.out_dir}/cts'
        self.ldsc_seg_dir = f'{self.out_dir}/ldsc_seg'
        self.start_ldsc_env=f'source /app/conda/bin/activate ldsc && cd {LDSC_bin_dir}'
        make_dir(self.out_dir)
        pass

    def make_gene_coor_by_gtf(self):
        def extract_gene_name(x):
            for k in str(x).split(';'):
                arr=re.split('\s+',k.strip())
                if arr[0]=='gene_name':
                    return arr[1].strip('"')
        df=pd.read_table(hg19_gtf,header=None,comment='#')
        df=df.loc[df[2]=='gene',]
        log(f'load {df.shape[0]} genes from {gene_annot_ldsc}')
        gene_name=df[8].map(extract_gene_name)
        codf=df[[0,3,4]]
        codf.columns=['CHR','START','END']
        codf['GENE']=gene_name
        codf.to_csv(gene_annot_ldsc,sep='\t',index=None)
        log(f'save gene coor to {gene_annot_ldsc}')

    def make_hmap3_snplist(self):
        df=pd.read_table(raw_snp_list_ldsc,header=0,dtype=str)
        df[['SNP']].to_csv(fine_hm3_snplist,index=None,header=None,lineterminator='\n')
        log(f'save to {fine_hm3_snplist}')

    def preprocess_tissue_specific_gene(self):
        wdir = pDESE_combine_rez_DIR
        top_gene_root_dir=self.top_gene_dir
        for f in os.listdir(wdir):
            if not f.endswith('.txt.gz'):
                continue
            log(f'{f} start')
            top_gene_dir= f'{top_gene_root_dir}/{f}'
            make_dir(top_gene_dir)
            raw_df = pd.read_table(f'{wdir}/{f}', index_col=0)
            df = extract_col(raw_df, '.z')
            with open(f'{top_gene_dir}/AllGene.txt', 'w') as bw:
                bw.write('\n'.join(df.index.tolist()))
            for column in df.columns:
                top_row_names = df[column].sort_values(ascending=False).head(LDSC_SEG_top_gene).index.tolist()
                with open(f'{top_gene_dir}/{column}.txt','w') as bw:
                    bw.write('\n'.join(top_row_names))
        pass


    def cal_tissue_ldscore(self):
        cmds = []
        for expr_name in os.listdir(self.top_gene_dir):
            edir=f'{self.top_gene_dir}/{expr_name}'
            for f in os.listdir(edir):
                # if re.search(r'AllGene|BrainCortex|ArteryAorta',f) is None:
                #     continue
                gene_set_path=f'{edir}/{f}'
                cell_ldsc_prefix=f'{self.tissue_ldsc_dir}/{expr_name}/{re.sub(".txt","",f)}'
                make_dir(os.path.dirname(cell_ldsc_prefix))
                for chr_i in range(1,23):
                    cmd=f'''
                         {self.start_ldsc_env} &&
                         python make_annot.py
                            --gene-set-file {gene_set_path}
                            --gene-coord-file {gene_annot_ldsc}
                            --windowsize 100000
                            --bimfile {eur_ref_plink}/1000G.EUR.QC.{chr_i}.bim
                            --annot-file {cell_ldsc_prefix}.{chr_i}.annot.gz
                        &&
                        python ldsc.py
                            --l2
                            --bfile {eur_ref_plink}/1000G.EUR.QC.{chr_i}
                            --ld-wind-cm 1
                            --annot {cell_ldsc_prefix}.{chr_i}.annot.gz
                            --thin-annot
                            --out {cell_ldsc_prefix}.{chr_i}
                            --print-snps {fine_hm3_snplist}
                    '''
                    cmd=re.sub('\s+',' ',cmd).strip()
                    cmds.append(cmd)
        batch_shell_plus(cmds,self.nt)

    def make_cts_file(self):
        make_dir(self.ldsc_cts_dir)
        for expr_name in os.listdir(self.tissue_ldsc_dir):
            edir=f'{self.tissue_ldsc_dir}/{expr_name}'
            expr_tag=f'{".".join(expr_name.split(".")[:2])}'
            cts_path=f'{self.ldsc_cts_dir}/{expr_tag}.ldcts'
            tissues = sorted(set(['.'.join(f.split('.annot.gz')[0].split('.')[:-1]) for f in os.listdir(edir) if
                                  f.endswith('.annot.gz')]))
            if 'AllGene' not in tissues:
                log(f'({expr_tag}) No background group, exit')
                continue
            with open(cts_path, 'w') as bw:
                for t in tissues:
                    if t=='AllGene':
                        continue
                    bw.write('\t'.join([t,f'{edir}/{t}.,{edir}/AllGene.'])+'\n')
            log(f'({expr_tag}) save {len(tissues)-1} tissues ldcts to {cts_path}')
        pass

    def make_gwas_raw(self):
        make_dir(LDSC_GWAS_dir)
        df=pd.read_excel(GWAS_raw_meta)
        cmds=[]
        for i in df.index:
            pheno=str(df.loc[i,"abbr"]).strip()
            # if pheno!='CAD':
            #     continue
            raw_path=f'{GWAS_raw_data_dir}/{str(df.loc[i,"category"]).strip()}/{str(df.loc[i,"dir"]).strip()}/{str(df.loc[i,"file_name"]).strip()}'
            if not os.path.isfile(raw_path):
                log(f'error: ({pheno}) no such file: {raw_path}')
                continue
            cmd = f'''
                 {self.start_ldsc_env} &&
                 python munge_sumstats.py
                    --chunksize 100000
                    --sumstats {raw_path} 
                    --merge-alleles {raw_snp_list_ldsc} 
                    --out {LDSC_GWAS_dir}/{pheno}.ldsc
                 '''
            cmd = re.sub('\s+', ' ', cmd).strip()
            snp_val=df.loc[i,'rsid']
            if not pd.isna(snp_val) and str(snp_val).strip()!='':
                cmd += f' --snp {str(snp_val).strip()}'
            n_sample=df.loc[i,'n_sample']
            if not pd.isna(n_sample) and str(n_sample).strip()!='':
                cmd += f' --N {int(n_sample)}'
            a1=df.loc[i,'a1']
            if not pd.isna(a1) and str(a1).strip()!='':
                cmd += f' --a1 {str(a1).strip()}'
            a2=df.loc[i,'a2']
            if not pd.isna(a2) and str(a2).strip()!='':
                cmd += f' --a2 {str(a2).strip()}'
            # log(f'({i}. {pheno}): {cmd}')
            cmds.append(cmd)
            # run_command(cmd)
            # return
        batch_shell_plus(cmds,self.nt)
        return

    def make_gwas(self):
        make_dir(LDSC_GWAS_dir)
        cmds=[]
        for f in os.listdir(Full_fine_GWAS_dir):
            raw_path=f'{Full_fine_GWAS_dir}/{f}'
            pheno=f.split('.')[0]
            cmd = f'''
                 {self.start_ldsc_env} &&
                 python munge_sumstats.py
                    --chunksize 100000
                    --sumstats {raw_path} 
                    --merge-alleles {raw_snp_list_ldsc} 
                    --out {LDSC_GWAS_dir}/{pheno}.ldsc
                    --N-col n_eff
                    --ignore beta,se
                 '''
            cmd = re.sub('\s+', ' ', cmd).strip()
            cmds.append(cmd)
            # run_command(cmd)
            # return
        batch_shell_plus(cmds,self.nt)
        return



    def run_regression(self):
        cmds=[]
        make_dir(self.ldsc_seg_dir)
        for gwas in os.listdir(LDSC_GWAS_dir):
            if not gwas.endswith('.sumstats.gz'):
                continue
            pheno=gwas.split('.')[0]
            for cts in os.listdir(self.ldsc_cts_dir):
                expr_tag=cts.split('.ldcts')[0]
                cmd = f'''
                     {self.start_ldsc_env} &&
                     python ldsc.py
                        --h2-cts {LDSC_GWAS_dir}/{gwas} 
                        --ref-ld-chr {LDSC_1kg_baseline} 
                        --out {self.ldsc_seg_dir}/{pheno}.{expr_tag} 
                        --ref-ld-chr-cts {self.ldsc_cts_dir}/{cts} 
                        --w-ld-chr {LDSC_weights_hm3} 
                     '''
                cmd = re.sub('\s+', ' ', cmd).strip()
                # log(f'({pheno} - {expr_tag}): {cmd}')
                # run_command(cmd)
                # return
                cmds.append(cmd)
        batch_shell_plus(cmds,self.nt)
def main():
    ldsc=LDSC(8)
    # ldsc.make_hmap3_snplist()
    # ldsc.make_gene_coor_by_gtf()
    # ldsc.preprocess_tissue_specific_gene()
    # ldsc.cal_tissue_ldscore()
    # ldsc.make_cts_file()
    # ldsc.make_gwas()
    # ldsc.run_regression()


if __name__ == '__main__':
    main()
