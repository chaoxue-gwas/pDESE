# -*- coding: utf-8 -*-
# @Author: Xue Chao
# @Time: 2024/07/01 11:32
# @Function: Analysis associated genes.

import urllib3
import os

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)
os.environ['HTTP_PROXY'] = 'socks5://127.0.0.1:7898'
os.environ['HTTPS_PROXY'] = 'socks5://127.0.0.1:7898'

import os
import pickle
import platform
import re
import shutil
import sys
import time

import numpy as np
import pandas as pd
import paramiko
import requests
import scipy
import seaborn
from expredix.projects.DriverNet.Result.AssociatedGene import __assocGeneVennPlot
from gprofiler import GProfiler
from matplotlib import pyplot as plt
import seaborn as sns
from venn import venn
import xml.dom.minidom as xmldom

from para import GWAS_ROOT, GWAS_case_dir, GWAS_meta_xlsx, output_DIR, COLORS_gene_compare_auc, \
    select_phenos, paper_out_dir, Pubmed_db_dir, Ecs_dir, pDESE_combine_rez_DIR, GWAS_raw_meta
from util import make_dir, read_line, log, get_gene_alias, replace_with_gene_symbol_in_refgene_of_df, \
    replace_with_gene_symbol_of_df, pvalue_adjust, run_command, LOCAL_DIR, batch_shell_plus


class NCBI:
    def __init__(self):
        self.sess=requests.session()
        pass
    def request_ncbi(self,url):
        res=self.sess.post(url)
        xobj=xmldom.parseString(res.text)
        count=xobj.documentElement.getElementsByTagName("Count")[0].firstChild.data
        pmids=[]
        for pid in xobj.documentElement.getElementsByTagName("Id"):
            pmids.append(str(pid.firstChild.data).strip())
        return count,pmids

    def single_trait_gene(self,traits,genes):
        geneTerm='+OR+'.join([f'({gene}[tiab])' for gene in genes])
        traitTerm='+OR+'.join([f'({t}[tiab])' for t in traits])
        base_url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?'
        query=f'db=pubmed&term=(({traitTerm})+AND+({geneTerm})+AND+(gene[tiab]+OR+genes[tiab]+OR+mRNA[tiab]+OR+protein[tiab]+OR+proteins[tiab]+OR+transcription[tiab]+OR+transcript[tiab]+OR+transcripts[tiab]+OR+expressed[tiab]+OR+expression[tiab]+OR+expressions[tiab]+OR+locus[tiab]+OR+loci[tiab]+OR+SNP[tiab]))&datetype=edat&retmax=100'
        url=f'{base_url}{query}'
        count=0
        pmids=[]
        while True:
            try:
                count,pmids=self.request_ncbi(url)
                break
            except:
                log(f'except! waiting for retrying ...')
                time.sleep(30)
                continue
        return str(count).strip(),pmids

    def batch_trait_gene(self,traits,genes,out_path):
        '''
        断点续传功能
        Gene同义符号搜索
        :param genes:
        :param out_path:
        :return:
        '''
        ## 获取基因同义词
        alias_genes=get_gene_alias(genes)
        if alias_genes is None:
            return None
        gene_alias_map={}
        for i in range(len(alias_genes)):
            gene_alias_map[genes[i]]=alias_genes[i]
        first = False
        egenes=[]
        if os.path.isfile(out_path):
            try:
                df=pd.read_table(out_path,header=0)
                egenes = list(df['gene'])
            except:
                egenes=[]
        else:
            first=True
        search_genes=[]
        for g in genes:
            if g not in egenes:
                search_genes.append(g)
        log(f'all: {len(genes)}, exist: {len(egenes)}, remain: {len(search_genes)}')
        with open(out_path,'a') as bw:
            if first:
                bw.write('\t'.join(['gene','count','pmids'])+'\n')
            i=0
            for g in search_genes:
                i+=1
                c,pids=self.single_trait_gene(traits,gene_alias_map[g])
                bw.write('\t'.join([g,c,','.join(pids)])+'\n')
                log(f'({i}/{len(search_genes)}) {g}: {c}')
                time.sleep(2)
class KGGSEE:
    def __init__(self, prefix):
        self.prefix = prefix
        self.kggsee_jar = f'{LOCAL_DIR}/lib/kggsee/kggsee.jar'
        self.resource_dir = f'{LOCAL_DIR}/lib/kggsee/resources'
        self.ref_eur_gty = f'{LOCAL_DIR}/lib/gty/EUR.hg19.vcf.gz'

    def gene_based_p_cutoff(self):
        p_cut = None
        i = 0
        for line in read_line(f'{self.prefix}.log'):
            if 'Significance level of p value cutoffs for the overall error rate' in line:
                i += 1
                p_cut = float(line.strip().split(':')[-1].strip())
        return p_cut

    def ecs_cond_assoc_genes(self,gene_score_name=''):
        Pcut = self.gene_based_p_cutoff()
        gp = f'{self.prefix}{gene_score_name}.gene.assoc.condi.txt'
        genes = []
        if os.path.isfile(gp):
            df = pd.read_table(gp)
            genes = list(df.loc[df['CondiECSP'] < Pcut, 'Gene'].values)
        return genes

    def ecs_assoc_genes(self,gene_score_name=''):
        gp = f'{self.prefix}{gene_score_name}.gene.assoc.condi.txt'
        genes = []
        if os.path.isfile(gp):
            df = pd.read_table(gp)
            genes = list(df['Gene'])
        return genes

    def ecs_cond_assoc_genes_df(self,gene_score_name=''):
        Pcut = self.gene_based_p_cutoff()
        gp = f'{self.prefix}{gene_score_name}.gene.assoc.condi.txt'
        df=None
        if os.path.isfile(gp):
            rdf = pd.read_table(gp)
            df = rdf.loc[rdf['CondiECSP'] < Pcut, :]
        return df

    def ecs_assoc_genes_df(self,gene_score_name=''):
        gp = f'{self.prefix}{gene_score_name}.gene.assoc.condi.txt'
        df=None
        if os.path.isfile(gp):
            df = pd.read_table(gp)
        return df


def getTPRandFPR(trueGenes,falseGenes,predTrueGenes,predFalseGenes):
    TP=len(set(predTrueGenes).intersection(set(trueGenes)))
    TN=len(set(predFalseGenes).intersection(set(falseGenes)))
    FP=len(set(predTrueGenes).intersection(set(falseGenes)))
    FN=len(set(predFalseGenes).intersection(set(trueGenes)))
    TPR,FPR=0,0
    if TP+FN!=0:
        TPR=TP/(TP+FN)
    if FP+TN!=0:
        FPR=FP/(FP+TN)
    return TPR,FPR,TP,len(predTrueGenes)

class PhenoGene:
    def __init__(self,db_name):
        self.db_name=db_name

    def access_db(self,phenotypes:[],out_tsv,genes):
        '''
        :return:
        '''
        support_dbs=['pubmed']
        if self.db_name not in support_dbs:
            raise Exception(f'Database: {self.db_name} is not supported! Supported list: {", ".join(support_dbs)}')
        make_dir(os.path.dirname(out_tsv))
        if self.db_name=='pubmed':
            NCBI().batch_trait_gene(phenotypes, genes, out_tsv)

    def eval_ROC_gene_based(self,sort_genes_dfs, ncbi_df, tags, libCutoff=5, title='', plot_step=1,print_info=False,plot_fig=True):
        colors=COLORS_gene_compare_auc
        tag_auc = {}
        if ncbi_df is None:
            return False, tag_auc
        db = ncbi_df
        if db.shape[0] == 0:
            return False, tag_auc
        trueGeneNum = db.loc[db['count'] >= libCutoff,].index.size
        rocVals = {}
        ppVal = {}
        k=-1
        for df in sort_genes_dfs:
            k+=1
            specDB = db.loc[[i for i in db.index if db.loc[i, 'gene'] in list(df['unique_gene'].values)],]
            trueGenes = list(specDB.loc[specDB['count'] >= libCutoff, 'gene'].values)
            falseGenes = [g for g in specDB['gene'] if g not in trueGenes]
            tag = tags[k]
            rocVals[tag] = []
            for i in np.arange(0,df.shape[0],plot_step):
                predTrueGenes = list(df.loc[df.index[:i], 'unique_gene'].values)
                predFalseGenes = list(df.loc[df.index[i:], 'unique_gene'].values)
                TPR, FPR, TP, PP = getTPRandFPR(trueGenes, falseGenes, predTrueGenes, predFalseGenes)
                rocVals[tag].append([TPR, FPR])
                # rocVals[tag].append([TP, PP])
            PP = df.index.size
            TP = len(set(df['unique_gene'].values).intersection(trueGenes))
            ppVal[tag] = [PP, TP]
        if print_info:
            print(f'all：{db.index.size}; positive：{trueGeneNum}；ratio：{trueGeneNum / db.index.size:.3f}')
        kc=-1
        xcor= {}
        ycor={}
        for k in tags:
            kc+=1
            xs = [z[1] for z in rocVals[k]]
            ys = [z[0] for z in rocVals[k]]
            auc = np.nansum(ys) / len(ys)
            if print_info:
                print(f'{k}: predict：{ppVal[k][0]}; TP：{ppVal[k][1]}；TPR：{ppVal[k][1] / ppVal[k][0]:.3f}')
            tag_auc[k] = auc
            xcor[k]=xs
            ycor[k]=ys
        # axe.plot([0, 1], [0, 1], '-', c='black')
        fig=None
        if plot_fig:
            fig, axe = plt.subplots(figsize=(3.5, 3.5))
            for kc in range(len(tags)):
                k=tags[kc]
                auc=tag_auc[k]
                axe.plot(xcor[k], ycor[k], label=f'{k}, AUC={auc:.3f}',c=colors[kc])
            axe.set_xlabel('FPR')
            axe.set_ylabel('TPR')
            axe.set_title(title)
            axe.spines['top'].set_visible(False)
            axe.spines['right'].set_visible(False)
            plt.legend(loc='lower right')
            # plt.show()
            plt.tight_layout()
        return fig, tag_auc


def compare_assoc_gene_AUC_plot(df, pheno_name):
    pale_cols=COLORS_gene_compare_auc
    data = {'pheno': [], 'auc': [], 'type': []}
    for c in df.columns:
        for i in df.index:
            data['pheno'].append(i)
            data['auc'].append(df.loc[i, c])
            data['type'].append(c)
    pldf = pd.DataFrame(data)
    if pldf.shape[0] == 0:
        return None
    units = df.shape[0]
    h_ratio = 6 / 4
    if units < 8:
        units = 8
        h_ratio = 4 / 2
    w = 6 / 12 * units
    h = w * h_ratio
    fig, ax = plt.subplots(figsize=(4, 3.3))
    sns.barplot(pldf, x='pheno', y='auc', hue='type', ax=ax,
                palette=pale_cols)  # palette=["#D77988", "#34ABBC"],
    ax.set_xticklabels(ax.get_xticklabels())  # rotation=70,ha='right'
    ax.set_ylim(0.5, np.max(df.values)+0.05)
    ax.set_xlabel('')
    ax.set_ylabel('AUC')
    ax.set_title(pheno_name)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    # plt.show()
    return fig

def compare_assoc_gene_AUC_group_plot(df,plot_box=False,two_side=True):
    test_h0='greater'
    if two_side:
        test_h0='two-sided'
    pales = COLORS_gene_compare_auc
    box_pale = pales[:df.shape[1]]
    data=[df[type].values for type in df.columns]
    types=df.columns.values
    fig, ax = plt.subplots(figsize=(1.15*len(data),3.4))
    if plot_box:
        seaborn.boxplot(data=data, ax=ax, palette=box_pale)
    else:
        seaborn.violinplot(data=data, ax=ax, palette=box_pale)
    ax.set_xticklabels(types)
    ## p value
    max_y = max([max(d) for d in data])
    min_y = min([min(d) for d in data])
    last_y=max_y*1.15
    y_high=max_y*1.35
    x_positions = ax.get_xticks()
    xi=x_positions[-1]
    for i in range(0 ,len(data)-1):
        i=len(data)-i-2
        xip1=x_positions[i]
        stat, pv = scipy.stats.ttest_rel(data[i], data[-1], alternative=test_h0)
        ax.plot((xip1,xi), (last_y, last_y), '-', c='black')
        ax.plot((xi, xi), (last_y, last_y*0.99), '-', c='black')
        ax.plot((xip1, xip1), (last_y, last_y * 0.99), '-', c='black')
        ax.text((xi + xip1) / 2, last_y*1.01, f'$P$ = {pv:.2g}', ha='center', va='bottom', fontsize=10, color='black')
        last_y+=max_y*0.09
    ax.set_ylim((min_y*0.7,y_high))
    ax.set_ylabel('AUC')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    return fig


def compare_assoc_gene_venn(genes,pheno_name, fig_path, show_size=True):
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    if show_size:
        fmt="{size}\n({percentage:.1f}%)"
    else:
        fmt="{percentage:.1f}%"
    venn(genes, ax=ax,fmt=fmt,cmap=COLORS_gene_compare_auc,legend_loc="lower right")
    plt.title(pheno_name,fontsize=20)
    plt.tight_layout()
    fig.savefig(fig_path)
    return True


def search_ncbi(ecs_dir,gwas_conf,db_local_dir):
    make_dir(db_local_dir)
    gdf = pd.read_excel(gwas_conf,dtype=str)
    gwas_alias = {}
    for i in gdf.index:
        gwas_alias[gdf.loc[i, 'abbr']] = f"{gdf.loc[i, 'ncbi_search_alias']}".split(';')
    phenos = sorted(set([x.split('.')[0] for x in os.listdir(ecs_dir) if os.path.isfile(f'{ecs_dir}/{x}')]))
    for pheno in phenos:
        out_path=f'{db_local_dir}/{pheno}.tsv'
        # if os.path.isfile(out_path):
        #     log(f'{pheno} exist, skip!')
        #     continue
        ## extract ecs genes
        genes = KGGSEE(f'{ecs_dir}/{pheno}').ecs_assoc_genes()
        genes=sorted(genes)
        log(f'start {pheno}: {gwas_alias[pheno]} with {len(genes)} genes')
        if len(genes)==0:
            continue
        pg = PhenoGene('pubmed')
        pg.access_db(gwas_alias[pheno],out_path,genes)

def eval_pheno_gene_ROC(kggsee_dir,ecs_dir,db_local_dir,fig_dir):
    make_dir(fig_dir)
    min_genes=50
    db_name='pubmed'
    phenos = select_phenos
    phe_auc={}
    # tags=[] # DGN
    tags=['P-value','RNA','Protein']
    kggsee_gene_scores=['CELL_2020.transcriptome.txt.gz','CELL_2020.proteome.txt.gz'] ## must match with tags. 'expr.mean.txt.gz',
    combine_data=[]
    combine_pheno_genes={}
    for pheno in phenos:
        log(f'start {pheno}')
        sig_genes = {}
        sort_genes_dfs = []
        genes_dfs=[]
        kgs=KGGSEE(f'{ecs_dir}/{pheno}')
        df=kgs.ecs_assoc_genes_df()
        sig_p_cutoff=kgs.gene_based_p_cutoff()
        genes_dfs.append([df,'CondiECSP','Gene'])
        if df is None or df.shape[0]<min_genes:
            log(f'ignore {pheno} due to genes < {min_genes}')
            continue
        for gene_score in kggsee_gene_scores:
            df = KGGSEE(f'{kggsee_dir}/{pheno}').ecs_assoc_genes_df(gene_score)
            genes_dfs.append([df,'CondiECSP','Gene'])
        for i in range(len(genes_dfs)):
            df, p_col, gene_col=genes_dfs[i]
            df['unique_p']=df[p_col]
            df['unique_gene']=df[gene_col].map(lambda x:f'{pheno}:{x}')
            df.sort_values(by='unique_p', inplace=True, ignore_index=True)
            sort_genes_dfs.append(df)
            sig_genes[tags[i]]=set(df.loc[df['unique_p']<sig_p_cutoff,'unique_gene'].values.tolist())

        ncbi_df=pd.read_table(f'{db_local_dir}/{pheno}.tsv')
        ncbi_df['gene']=ncbi_df['gene'].map(lambda x: f'{pheno}:{x}')

        combine_data.append([sort_genes_dfs,ncbi_df])

        fig, tag_auc = PhenoGene(db_name).eval_ROC_gene_based(sort_genes_dfs, ncbi_df,
                                          tags, 5,title=pheno)
        if fig:
            fig.savefig(f'{fig_dir}/{pheno}.{db_name}.fine.ROC.png')
            phe_auc[pheno] = [tag_auc[k] for k in tags]
        compare_assoc_gene_venn(sig_genes,pheno,f'{fig_dir}/{pheno}.fine.Venn.png')
        # plot go enrichment
        for k in tags:
            gs=[str(g).split(':')[1] for g in sig_genes[k]]
            go_enrich_plot(gs,f'{fig_dir}/{pheno}.{k}.go_enrich.png',title=f'{pheno} ({k})')
        for tag in sig_genes.keys():
            if tag not in combine_pheno_genes:
                combine_pheno_genes[tag]=[]
            combine_pheno_genes[tag]+=sorted([f'{pheno}:{sg}' for sg in sig_genes[tag]])
    combine_pheno_genes={k:[';'.join(combine_pheno_genes[k])] for k in combine_pheno_genes.keys()}
    pd.DataFrame(combine_pheno_genes).to_csv(f'{fig_dir}/compare.combine.pheno-gene.tsv',sep='\t',index=False)

    df = pd.DataFrame(phe_auc, index=tags).T
    fig = compare_assoc_gene_AUC_plot(df, '')
    fig.savefig(f'{fig_dir}/compare.AUC.{db_name}.fine.png')
    fig = compare_assoc_gene_AUC_group_plot(df)
    fig.savefig(f'{fig_dir}/compare.AUC.group.{db_name}.fine.box.png')
    df.to_csv(f'{fig_dir}/compare.AUC.{db_name}.fine.tsv',sep='\t')
    # combine ROC
    sdfss=[[] for i in range(len(tags))]
    ndfs=[]
    for i in range(len(combine_data)):
        sort_dfs,ncbi_df=combine_data[i]
        for k in range(len(sort_dfs)):
            sdfss[k].append(sort_dfs[k])
        ndfs.append(ncbi_df)
    com_sdfs=[pd.concat(dfs,ignore_index=True).sort_values(by=['unique_p']) for dfs in sdfss]
    com_ndf=pd.concat(ndfs,ignore_index=True)
    fig, tag_auc = PhenoGene(db_name).eval_ROC_gene_based(com_sdfs, com_ndf,
                                                          tags, 5, title="",plot_step=5)
    if fig:
        fig.savefig(f'{fig_dir}/combine.{db_name}.fine.ROC.png')
    with open(f'{fig_dir}/combine.{db_name}.fine.ROC.data.pyd','wb') as bw:
        pickle.dump((com_sdfs,com_ndf),bw)


def go_enrich_plot(m_genes, out_path,title='', save_max_term=10,plot_max_term=5, anno_dbs=['GO:BP', 'GO:CC', 'GO:MF']):
    out_tab=f'{out_path}.xlsx'
    if not os.path.isfile(out_tab):
        gp = GProfiler(return_dataframe=True)
        anno_df = gp.profile(organism='hsapiens', query=m_genes, sources=anno_dbs,
                             user_threshold=0.05)
        # anno_df=anno_df.loc[anno_df.sort_values(by=['p_value']).index[:max_term],:]
        remain_idx = []
        for cate, sdf in anno_df.groupby('source'):
            cut_term_n = save_max_term
            if cut_term_n > sdf.shape[0]:
                cut_term_n = sdf.shape[0]
            remain_idx += sdf.sort_values(by=['p_value']).index[:cut_term_n].tolist()
        anno_df = anno_df.loc[remain_idx,]
        anno_df = anno_df.sort_values(by=['source', 'p_value'])
        anno_df['p_value'] = -np.log10(anno_df['p_value'])
        anno_df.to_excel(out_tab)
        time.sleep(3)
    else:
        anno_df=pd.read_excel(out_tab)
    remain_idx = []
    for cate, sdf in anno_df.groupby('source'):
        cut_term_n = plot_max_term
        if cut_term_n > sdf.shape[0]:
            cut_term_n = sdf.shape[0]
        remain_idx += sdf.sort_values(by=['p_value'], ascending=False).index[:cut_term_n].tolist()
    anno_df = anno_df.loc[remain_idx,]
    try:
        bio_enrich_plot(anno_df, 'name', 'p_value', 'intersection_size', 'source',title=title, fig_path=out_path)
    except:
        log(f'WARNING: error in plot go enrichment')

def bio_enrich_plot(df:pd.DataFrame,y,x,size,color,title='',size_log=False,ax=None,fig_path=None,
                    max_letter=40,min_letter=40,adjusted_p=True):
    color_pale=['#E2614D','#FE9927','#4BAA52']
    # color_pale=['#06898A','#E47773']
    def process_string(s,min_l,max_l):
        if len(s) > max_l:
            s = s[:max_l] + '...'
        elif len(s) < min_l:
            s = ' ' * (min_l - len(s)) + s
        return s
    df=df.copy()

    df[y]=df[y].apply(lambda x:process_string(x,min_letter,max_letter))
    base_size=150
    size_legend_eg=[0.1,0.5,0.9]
    df['y_idx']=np.arange(df.shape[0],0,-1)
    raw_max_size=df[size].max()
    if size_log:
        df[size]=df[size].map(lambda x:np.log2(x+1))
    max_size=df[size].max()
    df['size_norm']=df[size]*base_size/max_size
    if ax is None:
        ratio=0.95
        h=df.shape[0]/16*3.7
        fig, ax = plt.subplots(figsize=(5,4))
    k=0
    for category, group in df.groupby(color):
        ax.scatter(group[x], group['y_idx'], s=group['size_norm'], label=f"{category}", alpha=1, c=color_pale[k])
        k+=1
    ax.set_yticks(df['y_idx'])
    ax.set_yticklabels(df[y])

    le1=ax.legend(loc='upper left',bbox_to_anchor=(1, 1))
    le1.set_title('Database')
    le1._legend_box.align = "left"
    le1.set_frame_on(False)
    for handle in le1.legendHandles:
        handle.set_sizes([base_size*0.7])
    p_text='P'
    if adjusted_p:
        p_text=r'Adjusted\ P'
    ax.set_xlabel(r'$-log_{10}('+p_text+r')$')
    ax2 = ax.twinx()
    for sl in size_legend_eg:
        show_eg=int(sl*max_size)
        eg=show_eg
        if size_log:
            show_eg=int(sl*raw_max_size//10*10)
            eg=np.log2(show_eg+1)
        ax2.scatter([], [], s=eg*base_size/max_size, label=f'{show_eg}', alpha=0.8, color='black')
    le2=ax2.legend(loc='upper left',bbox_to_anchor=(1, 0.64))
    le2.set_title('Count')
    le2._legend_box.align = "left"
    le2.set_frame_on(False)
    ax2.set_yticks([])
    ax.grid(axis='y', color='lightgrey', linestyle='--')
    ax.set_axisbelow(True)
    ax.grid(axis='x', color='lightgrey', linestyle='--')
    xmin, xmax = ax.get_xlim()
    new_xmin = xmin - 0.1 * (xmax - xmin)
    new_xmax = xmax + 0.1 * (xmax - xmin)
    ax.set_xlim(new_xmin, new_xmax)
    plt.title(f'{title}')
    plt.subplots_adjust(left=0.6, right=0.8, bottom=0.15, top=0.9)
    if fig_path is not None:
        make_dir(os.path.dirname(fig_path))
        plt.savefig(fig_path,dpi=200)
    plt.show()

def __gene_note(p_dgn,p_dese,p_cutoff,tags):
    tag_val=[f'{tags[0]} Sig. only',f'{tags[1]} Sig. only','Both Sig.','Unsig.']
    v=3
    if p_dgn<=p_cutoff and p_dese<=p_cutoff:
        v=2
    elif p_dgn<=p_cutoff and p_dese>p_cutoff:
        v=0
    elif p_dgn>p_cutoff and p_dese<=p_cutoff:
        v=1
    return v,tag_val[v]

def combine_assoc_genes(kggsee_dir,output_dir):
    phenos=select_phenos
    tags=['Protein','RNA']
    kggsee_gene_scores=['CELL_2020.proteome.txt.gz','CELL_2020.transcriptome.txt.gz']
    expr_tags=dict(zip(tags,kggsee_gene_scores))

    with pd.ExcelWriter(f'{output_dir}/associated_genes.xlsx') as bw:
        combine_pheno_genes = {}
        combine_tag_df = {}
        combine_ncbi_df = []
        combine_tag_auc = {}
        for pheno in phenos:
            print(pheno)
            ncbi_df = pd.read_table(f'{db_local_dir}/{pheno}.tsv')
            gene_pmids = {}
            for ni in ncbi_df.index:
                gene_pmids[ncbi_df.loc[ni, 'gene']] = ncbi_df.loc[ni, 'count']
            genes_dfs = []
            sig_genes = {}
            kp = f'{kggsee_dir}/{pheno}'
            ks = KGGSEE(kp)
            dfs = {}
            sig_p_cutoff = 1
            for m in tags:
                df = ks.ecs_assoc_genes_df(expr_tags[m])
                df.index=df['Gene']
                df = df.loc[~df.index.duplicated(), :]
                dfs[m] = df
                sig_p_cutoff = ks.gene_based_p_cutoff()
                genes_dfs.append([df, 'CondiECSP'])
            spec_cols = ['GeneScore', 'CondiECSP']
            com_df = dfs[f'{tags[0]}'].loc[:,
                     [x for x in dfs[f'{tags[0]}'].columns if x not in spec_cols]]
            sort_idx = com_df.index.values
            for m in tags:
                dfs[m] = dfs[m].loc[sort_idx, :]
            for scol in spec_cols:
                for m in tags:
                    com_df[f'{scol}({m})'] = dfs[m].loc[:, scol]
            note_tags = []
            note_v = []
            for g in com_df.index:
                tag_v, tag = __gene_note(com_df.loc[g, f'CondiECSP({tags[0]})'],
                                              com_df.loc[g, f'CondiECSP({tags[1]})'], sig_p_cutoff,tags)
                note_tags.append(tag)
                note_v.append(tag_v)
            com_df[f'note'] = note_tags
            com_df[f'note_val'] = note_v
            com_df['Pubmed count'] = com_df.index.map(lambda x: gene_pmids[x])
            com_df = com_df.sort_values(by=['note_val','Pubmed count', f'CondiECSP({tags[0]})'], ascending=[True, False, True])
            com_df = com_df.loc[:, [x for x in com_df.columns if x not in ['note_val','Gene']]]
            com_df.index.name = 'Gene'
            com_df.to_excel(bw, sheet_name=f'{pheno}')
            ## do
            db_name = 'pubmed'
            ncbi_df['gene'] = ncbi_df['gene'].map(lambda x: f'{pheno}:{x}')
            combine_ncbi_df.append(ncbi_df)
            sort_genes_dfs = []
            for i in range(len(tags)):
                tag = tags[i]
                df, p_col = genes_dfs[i]
                df['unique_p'] = df[p_col]
                df['unique_gene'] = df.index.map(lambda x: f'{pheno}:{x}')
                df.sort_values(by='unique_p', inplace=True, ignore_index=True)
                sort_genes_dfs.append(df)
                if tag not in combine_tag_df:
                    combine_tag_df[tag] = []
                combine_tag_df[tag].append(df)
                sig_genes[tag] = set(df.loc[df['unique_p'] < sig_p_cutoff, 'unique_gene'].values.tolist())
            fig, tag_auc = PhenoGene(db_name).eval_ROC_gene_based(sort_genes_dfs, ncbi_df,
                                                                  tags, 5, title=pheno, plot_fig=False)
            for i in range(len(tags)):
                tag = tags[i]
                if tag not in combine_tag_auc:
                    combine_tag_auc[tag] = []
                combine_tag_auc[tag].append(tag_auc[tag])
            for tag in sig_genes.keys():
                if tag not in combine_pheno_genes:
                    combine_pheno_genes[tag] = []
                combine_pheno_genes[tag] += sorted([f'{pheno}:{sg}' for sg in sig_genes[tag]])
        ## plot auc volint
        auc_df = pd.DataFrame(combine_tag_auc)
        fig = compare_assoc_gene_AUC_group_plot(auc_df, two_side=True)
        fig.savefig(f'{output_dir}/compare_gene_AUC_group.png')
        plt.close()
        ddfs = []
        for m in tags:
            ddf = pd.concat(combine_tag_df[m], axis=0, ignore_index=True)
            ddf = ddf.sort_values(by=['unique_p'])
            ddfs.append(ddf)
        fn_df = pd.concat(combine_ncbi_df, axis=0, ignore_index=True)
        fig, tag_auc = PhenoGene(db_name).eval_ROC_gene_based(ddfs, fn_df,
                                                              tags, 5, title='', plot_fig=True)
        fig.savefig(f'{output_dir}/compare_gene_ROC.png')


def check_expr(genes,tissues):
    tags=['Protein','RNA']
    kggsee_gene_scores=['CELL_2020.proteome.REZ.webapp.genes.txt.gz','CELL_2020.transcriptome.REZ.webapp.genes.txt.gz']
    expr_tags=dict(zip(tags,kggsee_gene_scores))
    wdir=pDESE_combine_rez_DIR
    dfs={}
    for tag in tags:
        df=pd.read_table(f'{wdir}/{expr_tags[tag]}',index_col=0)
        dfs[tag]=df
    score={}
    clean_genes=[]
    for g in genes:
        has_g=True
        for tag in tags:
            if g not in dfs[tag].index:
                has_g=False
                break
        if has_g:
            clean_genes.append(g)

    for gene in clean_genes:
        score[gene]={}
        for t in tissues:
            tagv={}
            for tag in tags:
                df = dfs[tag]
                vals = {}
                cols=[c for c in df.columns if str(c).split('.')[0]==t]
                v=df.loc[gene,cols]
                for id in v.index:
                    ka=str(id).split('.')[-1]
                    va=v[id]
                    if str(ka).__contains__('rank*'):
                        fm=int(str(ka).strip().split('*')[1])
                        all=12230
                        fz=va
                        fz=fz+all-fm
                        ka='rank'
                        va=fz/all
                    vals[f'{ka}']=va
                tagv[tag]=vals
            sv={}
            for k in vals.keys():
                sv[k]=[]
                for tag in tags:
                    sv[k].append(tagv[tag][k])
            score[gene][t]=sv
    return score

def stat_protein_prio_gene(output_dir):
    pheno_tissues={'BIP':['BrainCortex','BrainCerebellum'],'SCZ':['BrainCortex','BrainCerebellum'],
                   'CAD':['ArteryCoronary'],'CD':['Spleen'],'RA':['Spleen']}
    result_excel=f'{output_dir}/associated_genes.xlsx'
    dfs=pd.ExcelFile(result_excel)
    new_cols=['Disease','Gene','P(Protein)','P(RNA)','Pubmed count','Associated tissue','Rank(Protein)','Rank(RNA)','Diff. rank']
    data=[]
    re_cols=['CondiECSP(Protein)','CondiECSP(RNA)','Pubmed count']
    for pheno in dfs.sheet_names:
        tissues=pheno_tissues[pheno]
        df=dfs.parse(pheno,index_col=0)
        gs=df.loc[df['note']=='Protein Sig. only',:].index.to_list()
        score=check_expr(gs,tissues)
        # print(score)
        for g in score.keys():
            common_vals=[]
            for c in re_cols:
                common_vals.append(df.loc[g,c])
            for t in score[g].keys():
                vals = [pheno, g]+common_vals+[t]
                if 'rank' not in score[g][t]:
                    continue
                rank_arr=score[g][t]['rank']
                vals+=rank_arr
                vals.append(rank_arr[0]-rank_arr[1])
                data.append(vals)
    fdf=pd.DataFrame(data,columns=new_cols)
    fdf.to_excel(f'{output_dir}/assoc_gene.analysis.xlsx',index=False)

if __name__ == '__main__':
    kggsee_dir = f'{output_DIR}/dese'
    db_local_dir=Pubmed_db_dir
    ecs_dir=Ecs_dir
    fig_dir=f'{paper_out_dir}/fig_gene'
    # search in NCBI PubMed.
    search_ncbi(ecs_dir,GWAS_raw_meta,db_local_dir)
    # evaluate
    eval_pheno_gene_ROC(kggsee_dir,ecs_dir,db_local_dir,fig_dir)
    # combine
    combine_assoc_genes(kggsee_dir,fig_dir)
    # analyze disease-gene associations only identified by protein-based approach
    stat_protein_prio_gene(fig_dir)