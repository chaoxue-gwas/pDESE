# -*- coding: utf-8 -*-
# @Author: Xue Chao
# @Time: 2024/07/01 11:32
# @Function: Analysis associated tissues.
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn

from paper.para import output_DIR, paper_out_dir
from paper.util import make_dir

ProteomeDBs=['CELL_2020.proteome','CELL_2020.transcriptome']
db_show={'CELL_2020.proteome':'protein','CELL_2020.transcriptome':'RNA'}
ExcludeGWAS=['PTSD','TS','AN','SMK']
Methods=['ldsc','magma','dese']
make_dir(paper_out_dir)


def common_combine_assoc_tissue_results(results_dir,dbs,method,log_suffix,assoc_tissue_suffix,p_col,no_point=False,
                                        comment=None,sep='\s+',skipfooter=0):
    phenos=sorted(set([f.split('.')[0] for f in os.listdir(results_dir) if f.endswith(log_suffix)]))
    phenos=[p for p in phenos if p not in ExcludeGWAS]
    for db in dbs:
        dfs=[]
        for p in phenos:
            tpath=f'{results_dir}/{p}.{db}.{assoc_tissue_suffix}'
            if no_point:
                tpath=f'{results_dir}/{p}{db}.{assoc_tissue_suffix}'
            df=pd.read_table(tpath,sep=sep,index_col=0,comment=comment,skipfooter=skipfooter,engine='python')
            sdf=df[[p_col]]
            sdf.columns=[p]
            dfs.append(sdf)
        fdf=pd.concat(dfs,axis=1)
        fdf.sort_index(inplace=True)
        fdf=fdf.applymap(lambda x:0 if pd.isna(x) else -np.log10(x))
        fdf.index.name='Tissue'
        fdf.to_excel(f'{paper_out_dir}/{method}.{db}.assoc_tissue.xlsx')

def combine_results():
    ldsc_dir = f'{output_DIR}/ldsc/ldsc_seg'
    common_combine_assoc_tissue_results(ldsc_dir,ProteomeDBs,'ldsc','.log','cell_type_results.txt','Coefficient_P_value')
    magma_dir = f'{output_DIR}/magma'
    common_combine_assoc_tissue_results(magma_dir, ProteomeDBs, 'magma', '.log', 'gsa.out',
                                        'P',comment='#')
    dese_dir = f'{output_DIR}/dese'
    common_combine_assoc_tissue_results(dese_dir, ProteomeDBs, 'dese', '.log', 'txt.gz.celltype.txt',
                                        'Unadjusted(p)',no_point=True,sep='\t',skipfooter=1)


def global_show():
    for m in Methods:
        for db in ProteomeDBs:
            df=pd.read_excel(f'{paper_out_dir}/{m}.{db}.assoc_tissue.xlsx',index_col=0)
            g = seaborn.clustermap(df,cmap="Reds",vmax=5,vmin=0, figsize=(9, 9))
            plt.show()

def bar_show():
    method_lab={'ldsc':'S-LDSC','magma':'MAGMA','dese':'DESE'}
    fig_dir=f'{paper_out_dir}/fig_fine'
    make_dir(fig_dir)
    results_dir=f'{output_DIR}/ldsc/ldsc_seg'
    phenos=sorted(set([f.split('.')[0] for f in os.listdir(results_dir) if f.endswith('.log')]))
    select_phenos=[p for p in phenos if p not in ExcludeGWAS]
    select_phenos=['MDD','BIP','SCZ','CAD','AF','CD','MS','RA','SLE','T1D']

    for phe in select_phenos:
        for db in ProteomeDBs:
            dfs=[]
            for m in Methods:
                df=pd.read_excel(f'{paper_out_dir}/{m}.{db}.assoc_tissue.xlsx',index_col=0)
                df=df[[phe]]
                df.columns=[m]
                dfs.append(df)
            fdf=pd.concat(dfs,axis=1)
            fdf['CombineRank']=fdf.rank(method="average", ascending=True).mean(axis=1)
            fdf['Tissue']=fdf.index
            fdf = fdf.sort_values("CombineRank", ascending=False)

            seaborn.set_theme(style="whitegrid")
            seaborn.set_theme(style="white")
            fig, ax1 = plt.subplots(figsize=(8, 4))
            ax1 = seaborn.barplot(x="Tissue", y="CombineRank", data=fdf, color="skyblue", edgecolor="black")
            ax1.set_xlabel("")
            ax1.set_ylabel("Mean of Ranks", color="black")
            ax1.tick_params(axis='y', labelcolor="black")
            plt.xticks(rotation=90)
            ax2 = ax1.twinx()
            ax2.set_ylabel("-log10(P)", color="black")

            markers = ['o','s','D']
            marker_colors = ['#F47381','#D58C2A','#53AC34']

            for i, method in enumerate(['ldsc', 'magma', 'dese']):
                ax2.scatter(fdf["Tissue"], fdf[method], label=method_lab[method], marker=markers[i], color=marker_colors[i],
                            s=30)
            ax2.axhline(-np.log10(0.05/fdf.shape[0]), color='gray', linestyle='dashed', linewidth=1)
            ax2.legend(title="Methods", loc="upper right")
            ax2.tick_params(axis='y', labelcolor="black")

            plt.title(f"{phe} ({db_show[db]})")
            plt.tight_layout()
            # plt.show()
            plt.savefig(f'{fig_dir}/{phe}.{db}.assoc_tissue.jpg',dpi=300)
        # return


def main():
    combine_results()
    global_show()
    bar_show()

if __name__ == '__main__':
    main()