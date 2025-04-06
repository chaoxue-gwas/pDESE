# -*- coding: utf-8 -*-
# @Author: Xue Chao
# @Time: 2024/05/24 16:09
# @Function:
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats
import seaborn
from matplotlib.patches import Patch
from scipy import stats

from paper.para import pDESE_combine_data_DIR, pDESE_combine_rez_DIR, pDESE_combine_meta, paper_out_dir
from paper.util import extract_col, cluster_df, generate_colormap, make_dir


def load_tissue_meta():
    df=pd.read_excel(pDESE_combine_meta)
    df=df.applymap(lambda x:str(x).strip())
    meta_dict={'raw':{},'cate':{}}
    for i in df.index:
        meta_dict['raw'][df.loc[i,'formatted']]=df.loc[i,'raw']
        meta_dict['cate'][df.loc[i, 'raw']] = df.loc[i, 'category_abbr']
    return meta_dict



def proteome_corr():
    wdir=pDESE_combine_rez_DIR
    meta = load_tissue_meta()
    for f in os.listdir(wdir):
        if not f.endswith('.txt.gz'):
            continue
        print(f)
        raw_df=pd.read_table(f'{wdir}/{f}',index_col=0)
        df=extract_col(raw_df,'.z')
        df.columns=df.columns.map(lambda x:meta['raw'][x])
        cc_df=df.corr(method='spearman')
        make_dir(paper_out_dir)
        cc_df.to_excel(f'{paper_out_dir}/protome.corr.xlsx')
        # cc_df=cluster_df(cc_df)
        # sorted_index=sorted(cc_df.index.to_list(),key=lambda x:meta[x][1])
        # cc_df=cc_df.loc[sorted_index,:].loc[:,sorted_index]
        cates=cc_df.index.map(lambda x:meta['cate'][x])
        row_colors,row_palette =generate_colormap(cates)
        col_colors=row_colors
        g = seaborn.clustermap(cc_df,cmap="bwr",center=0,vmax=1,vmin=-1,
                           row_colors=row_colors, col_colors=col_colors, figsize=(9, 9))
        colorbar = g.ax_heatmap.collections[0].colorbar
        colorbar.set_label("Spearman R")
        fig_legend, ax_legend = plt.subplots(figsize=(3, 3))
        row_legend_patches = sorted([Patch(color=color, label=label) for label, color in row_palette.items()],
                                    key=lambda x: x.get_label())
        ax_legend.set_frame_on(False)
        ax_legend.legend(handles=row_legend_patches, title="Tissue Categories", loc="center")
        ax_legend.set_xticks([])
        ax_legend.set_yticks([])
        ax_legend.axis("off")
        plt.show()
        plt.close()

        return

def proteome_trans_corr_hist():
    wdir=pDESE_combine_rez_DIR
    prefixs=['CELL_2020.proteome','CELL_2020.transcriptome']
    raw_dfs=[]
    meta = load_tissue_meta()
    for k in prefixs:
        raw_df=pd.read_table(f'{wdir}/{k}.REZ.webapp.genes.txt.gz',index_col=0)
        df=extract_col(raw_df,'.z')
        df.columns=df.columns.map(lambda x:meta['raw'][x])
        raw_dfs.append(df)
    tissues=raw_dfs[0].columns
    genes=raw_dfs[0].index
    dfs=[]
    for df in raw_dfs:
        dfs.append(df[tissues].loc[genes,:])
    ccs=[]
    ps=[]
    for t in tissues:
        cc,p=scipy.stats.spearmanr(dfs[0][t],dfs[1][t],nan_policy='omit')
        ccs.append(cc)
        ps.append(p)
    cates=[meta['cate'][t] for t in tissues]
    ylab='Spearman R'
    hdf=pd.DataFrame({'Tissue':tissues,'Category':cates,ylab:ccs,'P':ps})
    plt.figure(figsize=(6, 8))
    hdf = hdf.sort_values(by=["Category", "Tissue"], ascending=[True, True], ignore_index=True)
    tc, type_colors = generate_colormap(hdf['Category'])
    ax = seaborn.barplot(y='Tissue', x=ylab, data=hdf, hue='Category', palette=type_colors)
    mean_value = hdf[ylab].mean()
    confidence_interval = stats.t.interval(0.95, len(hdf[ylab]) - 1, loc=mean_value, scale=hdf[ylab].sem())
    plt.axvline(mean_value, color='gray', linestyle='dashed', linewidth=1)
    # plt.text(mean_value + 1, -1, f'Mean: {mean_value:.2f}', color='gray', verticalalignment='bottom')
    # for index, row in hdf.iterrows():
    #     if row['P'] < 0.01:
    #         ax.text(row[ylab], index, '*', color='red', ha='left', va='center', fontsize=12, fontweight='bold')
    plt.legend(title="Category", bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.xlabel(ylab)
    x_min, x_max = plt.xlim()
    plt.xlim(x_min, x_max * 1.2)
    y_min, y_max = plt.ylim()
    plt.ylim(y_min+0.3, y_max - 0.3)
    # plt.ylabel("Cell")
    plt.title("")
    plt.subplots_adjust(left=0.3, right=0.7, top=0.9, bottom=0.1)
    plt.show()
    print(f'Mean: {mean_value:.4f}; 95% CI: {confidence_interval}')
    hdf.to_excel(f'{paper_out_dir}/protome_vs_transcriptome.corr.xlsx',index=False)


def main():
    proteome_corr()
    proteome_trans_corr_hist()


if __name__ == '__main__':
    main()
