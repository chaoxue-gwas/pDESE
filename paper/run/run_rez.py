# -*- coding: utf-8 -*-
# @Author: Xue Chao
# @Time: 2025/02/27 14:56
# @Function: Calculate tissue-specific gene expression by REZ (robust z score, Genome Biology, 2019).
import os

from paper.para import pDESE_combine_data_DIR, pDESE_combine_rez_DIR
from paper.util import kggsee_rez, LOCAL_DIR, run_command


def cal_rez():
    kggsee_dir = f'{LOCAL_DIR}/lib/kggsee'
    kggsee_jar = f'{kggsee_dir}/kggsee.jar'
    resource_dir = f'{kggsee_dir}/resources'
    gene_score_dir=pDESE_combine_data_DIR
    gene_rez_dir=pDESE_combine_rez_DIR
    for f in os.listdir(gene_score_dir):
        expr_path=f'{gene_score_dir}/{f}'
        out_path=f'{gene_rez_dir}/{f.split(".txt.gz")[0]}'
        cmd=kggsee_rez(expr_path,out_path,5,kggsee_jar,resource_dir)
        run_command(cmd)
    pass

if __name__ == '__main__':
    cal_rez()