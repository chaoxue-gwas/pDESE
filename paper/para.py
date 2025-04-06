# -*- coding: utf-8 -*-
import platform

# Common
DATA_DIR = '/home/xc/local/data'
KGGSEE_RESOURCE = '/home/xc/local/data/resources/ToolResource/kggsee/resources'
if platform.system() == 'Windows':
    DATA_DIR = r'E:\WorkData\syncHPC\home\data'
    KGGSEE_RESOURCE=r'E:\WorkData\syncHPC\home\data\resources\ToolResource\kggsee\resources'

# Tool resource
LDSC_resource_dir=f'{DATA_DIR}/resources/ToolResource/LDSC'

# Proteome
PROTEOME_DIR = f'{DATA_DIR}/resources/Proteome'
RAW_HPA_PATH = f'{PROTEOME_DIR}/HPA/raw/normal_tissue.tsv'
RAW_CELL_2020_PRO_PATH = f'{PROTEOME_DIR}/Cell_2020/raw/cell_2020.proteome.xlsx' ## sheet: 'E protein tissue median'
RAW_CELL_2020_TRA_PATH = f'{PROTEOME_DIR}/Cell_2020/raw/cell_2020.transcriptome.xlsx' ## sheet: 'B RNA tissue median'
CELL_2020_PRO_PATH = f'{PROTEOME_DIR}/Cell_2020/analysis/gene_score/proteome.normalized.median.txt.gz'
CELL_2020_TRA_PATH = f'{PROTEOME_DIR}/Cell_2020/analysis/gene_score/transcriptome.log2_tpm.median.txt.gz'
# CELL_2020_PRO_NORM_PATH = f'{PROTEOME_DIR}/Cell_2020/analysis/norm_gene_score/proteome.txt.gz'
# CELL_2020_TRA_NORM_PATH = f'{PROTEOME_DIR}/Cell_2020/analysis/norm_gene_score/transcriptome.txt.gz'

pDESE_combine_data_DIR=f'{PROTEOME_DIR}/pDESE_combine/norm_gene_score'
pDESE_combine_rez_DIR=f'{PROTEOME_DIR}/pDESE_combine/norm_gene_rez'
pDESE_combine_data_noNA_DIR=f'{PROTEOME_DIR}/pDESE_combine/norm_gene_score_no_na'
pDESE_combine_data_magma_DIR=f'{PROTEOME_DIR}/pDESE_combine/norm_magma'

pDESE_combine_meta=f'{PROTEOME_DIR}/pDESE_combine/meta.xlsx'

## super dir

# Parameters for analysis in paper.
META_DIR=f'{DATA_DIR}/resources/meta'
transcriptome_meta_xlsx=f'{META_DIR}/transcriptome.xlsx'
GWAS_meta_xlsx=f'{META_DIR}/GWAS.xlsx'
meta_transcriptome_dir=f'{META_DIR}/transcriptome'

# GWAS
GWAS_ROOT=f'{DATA_DIR}/resources/GWAS'
GWAS_case_dir=f'{GWAS_ROOT}/gold_case/formatted'
GWAS_raw_data_dir=f'{GWAS_ROOT}/gold_case/raw'
GWAS_raw_meta=f'{GWAS_raw_data_dir}/meta_v2.xlsx'
LDSC_GWAS_dir=f'{GWAS_ROOT}/gold_case/ldsc_formatted'
Full_fine_GWAS_dir=f'{GWAS_ROOT}/gold_case/pmglab_formatted'
MAGMA_GWAS_dir=f'{GWAS_ROOT}/gold_case/magma_formatted'

# SNP annotation, rsid
SNP_rsid_dir=f'{DATA_DIR}/resources/VariantAnnotation/rsid'
raw_hg19_rsid_path=f'{SNP_rsid_dir}/human_9606_b151_GRCh37p13.common_all_20180423.vcf.gz'
raw_hg38_rsid_path=f'{SNP_rsid_dir}/human_9606_b151_GRCh38p7.common_all_20180418.vcf.gz'
hg19_rsid_path=f'{SNP_rsid_dir}/hg19.common.vcf.gz'
hg38_rsid_path=f'{SNP_rsid_dir}/hg38.common.vcf.gz'

# Output
output_main_version='pdese_0227'
output_relative_dir=f'projects/{output_main_version}'
output_DIR=f'{DATA_DIR}/{output_relative_dir}' #
## Paper
paper_out_dir=f'{output_DIR}/paper'

# Gene annotation
## from https://github.com/broadinstitute/gtex-pipeline/blob/master/TOPMed_RNAseq_pipeline.md
gtexv8_pip_gene_annot_gtf=f'{DATA_DIR}/resources/GeneAnnotation/gencode.v34.annotation.gtf.gz'
hg19_gtf=f'{DATA_DIR}/resources/GeneAnnotation/Homo_sapiens.GRCh37.87.gtf.gz'
## ENSG annotation for LDSC-SEG
gene_annot_ldsc=f'{LDSC_resource_dir}/GRCh37.ENSG_coord.txt'
raw_snp_list_ldsc=f'{LDSC_resource_dir}/w_hm3.snplist.gz'
fine_hm3_snplist=f'{LDSC_resource_dir}/fine.w_hm3.snplist.gz'
eur_ref_plink=f'{LDSC_resource_dir}/1000G_Phase3_plinkfiles'
LDSC_1kg_baseline=f'{LDSC_resource_dir}/1000G_Phase3_baseline_v1.2/baseline.'
LDSC_weights_hm3=f'{LDSC_resource_dir}/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.'

## MAGMA resources
magma_resource_dir=f'{DATA_DIR}/resources/ToolResource/MAGMA'
magma_gene_annot = f'{magma_resource_dir}/NCBI37.3/NCBI37.3.gene.loc'
magma_ref_geno = f'{magma_resource_dir}/g1000_eur/g1000_eur'

COLORS_gene_compare_auc=['#00bf7d','#FFA404','#F8766D'] # #00BFC4
CELL_CATE_colors= ['#DB5F57','#56DB5F','#5F56DB']

select_phenos = ['BIP', 'SCZ', 'CAD', 'CD', 'RA']

## gene based resource
Pubmed_db_dir=f'{GWAS_ROOT}/gold_case/gene_based/pubmed'
Ecs_dir=f'{GWAS_ROOT}/gold_case/gene_based/ecs'

