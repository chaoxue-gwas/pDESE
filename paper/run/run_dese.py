# -*- coding: utf-8 -*-
# @Author: Xue Chao
# @Time: 2025/03/12 10:16
# @Function: run pipeline to estimate associated tissue and genes by DESE.

import os
import re
import sys

from paper.para import Full_fine_GWAS_dir, output_DIR, pDESE_combine_data_DIR, GWAS_case_dir, \
    pDESE_combine_data_noNA_DIR
from paper.util import kggsee_dese, LOCAL_DIR, batch_shell, make_dir


def DESE():
    remove_hlc=True
    PSC_ROOT='/webapp/psc'
    java = '/usr/local/bin/java8'
    kggsum_dir = f'{PSC_ROOT}/kggsum_web'
    kggsum_jar = f'{kggsum_dir}/kggsum20250312.jar'
    ref_genotype_dir = f'{kggsum_dir}/psc_resources/ref_genotype/1KG_super_pop_combine'
    ntasks=int(sys.argv[1])
    gwas_dir=GWAS_case_dir
    gene_score_dir=pDESE_combine_data_noNA_DIR
    out_dir=f'{output_DIR}/dese'
    make_dir(out_dir)
    score_files = sorted([f for f in os.listdir(gene_score_dir)])
    gene_score_files = ' '.join([f'--gene-score-file file={gene_score_dir}/{sf} calcSpecificity=y' for sf in score_files])
    cmds = []
    for phe in os.listdir(gwas_dir):
        phe_abbr=phe.split('.')[0]
        phe_gwas_path=f'{gwas_dir}/{phe}'
        ref_ver = 'hg19'
        ref_vcf_ver = 'hg38'
        pop = 'EUR'
        mhc_size = 'chr6:28477797~33448354'
        if ref_ver == 'hg38':
            mhc_size = 'chr6:28510120~33480577'
        mhc_para = ''
        if remove_hlc:
            mhc_para = f'exclude={mhc_size}'
        custom_para = f'''
            assoc
            --ref-gty-file {ref_genotype_dir}/{pop}.{ref_vcf_ver}.vcf.bgz
             refG={ref_vcf_ver}
            --sum-file {phe_gwas_path}
             cp12Cols=chr,bp
             pbsCols=p
             refG={ref_ver} {mhc_para}
            --threads 8
            {gene_score_files}
            --output {out_dir}/{phe_abbr}
            --gene-p-cut 0.05
            --gene-multiple-testing bhfdr
            --max-condi-gene 1000
        '''
        ### run DESE
        clean_para = re.sub('\s+', ' ', custom_para)
        cmd = f'cd {kggsum_dir} && {java} -Xmx64G -jar {kggsum_jar} {clean_para}'
        cmds.append(cmd)
        break
    batch_shell(cmds,ntasks,f'{out_dir}/batch_run.log')

if __name__ == '__main__':
    DESE()
