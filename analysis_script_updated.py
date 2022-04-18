'''
Distal Fetal Lung CELL ANALYSIS SCRIPT
By Josh Wu
13 April, 2020


Contains of Distal Lung for Renee Conway
'''
import sys
sys.path.insert(0,'C:/Users/Josh/Desktop/sca_run/api/')
#from scanpy_spence import *
from classes.sca_params import *
import SCARunner as SCARunner

figdir = './figures_062921_dpt/'

an_params = sca_params()

an_params.storage_mount_point = 'Z:/'

## IDs of samples as represented in the metadata table
# an_params.sample_list = ['2268-2','2598-33','2321-2','2288-4','2444-1'] # Distal Fetal Lung data

an_params.sample_list = ['2288-4']

an_params.add_gene_list(markers = ['EGR1','RSPO2','FGFR4','WNT2','CDO1','BMP5','LIFR',
								   'COL1A1','COL1A2','ELN','SPARC','MGP','FBLN1','SOX11',
								   'PDGFRA','PTCH1','HOXB5','HOXB4','NGFR','FGF7','BMP4',
								   'TAGLN','ACTA2','MYL9','MYH11','ACTG2','HHIP','FOXF1',
								   'WNT5A'], 
						label='mesenchyme_list')

an_params.add_gene_list(markers = ['VIM','POSTN','DCN','TCF21','COL1A2','COL3A1','MKI67',
								   'TOP2A','CDK1','TYMS','H2AFZ','STMN1','EPCAM','KRT18',
								   'KRT8','CLDN6','TAGLN','ACTA2','PDGFRB','MGP','CDH5',
								   'CLDN5','ESAM','KDR','FLT1','PTPRC','CD37','CORO1A',
								   'LCP1','CD53','LAPTM5','S100B','ELAVL4','TUBB2B','STMN2',
								   'ASCL1','NNAT','GRP','MPZ'],
					 label='general_cell_type_list')


an_params.add_gene_list(markers = ['RGS13','RPS4Y1','CBLN2','APOC1','RPS4Y1','SOD3','CA3',
								   'ANXA1','SDC2','STC1','SYT1','IGSF1','KCNK17','NKX3-2',
								   'NTRK2','CHODL','CRYM','GATA2','PRSS23','SERPINF1'],
						label='ranked_mesenchyme_list')

an_params.add_gene_list(markers = ['RPS4Y1','RSPO2','WNT2','CDO1','BMP5','LIFR','FGFR4','EGR1',
								  'GUCA1A','CA3','SOX11','PDGFRA','PTCH1','HOXB4','SDC2','STC1',
								  'MAPK14','KCNK17','ELN','MGP','BGN','PRSS23','SERPINF1','CRYM',
								  'AGTR2','FGF7','BMP4','TAGLN','ACTA2','MYL9','MYH11','ACTG2',
								  'HHIP','FOXF1','WNT5A'],
					   label='updated_mesenchyme_3')

an_params.add_gene_list(markers = ['EPCAM','KRT18','KRT8','CLDN6','VIM','POSTN','DCN','TCF21',
								   'COL1A2','COL3A1','TAGLN','ACTA2','PDGFRB','RGS5','PTPRC',
								   'CD37','CORO1A','LCP1','CD53','LAPTM5','CDH5','CLDN5',
								   'ESAM','KDR','FLT1','S100B','ELAVL4','TUBB2B','STMN2',
								   'ASCL1','NNAT','GRP','MPZ','MKI67','TOP2A','CDK1','TYMS',
								   'H2AFZ','STMN1','LGR5'],
						label='general_cell_type_list_updated')

an_params.add_gene_list(markers = ['VIM','POSTN','DCN','TCF21','COL1A2','COL3A1','MKI67',
								   'TOP2A','CDK1','TYMS','H2AFZ','STMN1','EPCAM','KRT18',
								   'KRT8','CLDN6','TAGLN','ACTA2','PDGFRB','MGP','CDH5',
								   'CLDN5','ESAM','KDR','FLT1','PTPRC','CD37','CORO1A',
								   'LCP1','CD53','LAPTM5','S100B','ELAVL4','TUBB2B',
								   'STMN2','ASCL1','NNAT','GRP','MPZ','CHGA','SYN1','LGR5'],
						label='general_cell_type_list_explants')

an_params.add_gene_list(markers = ['WNT1','WNT2','WNT2B','WNT3','WNT3A','WNT4','WNT5A',
								   'WNT5B','WNT6','WNT7A','WNT7B','WNT8A','WNT8B','WNT9A',
								   'WNT9B','WNT10A','WNT10B','WNT11','WNT16','DVL1','DVL2',
								   'DVL3','HNF1A','TCF3','TCF4','LEF1','LRP5','LRP6',
								   'FZD1','FZD2','FZD3','FZD4','FZD5','FZD6','FZD7','FZD8',
								   'FZD9','FZD10','ANAPC1','APC2','DKK1','DKK2','DKK3',
								   'DKK4','GSK3A','GSK3B','PORCN','RSPO1','RSPO2','RSPO3',
								   'RSPO4','LGR4','LGR5','LGR6','RNF43','ZNRF3','AXIN2'],
					 	label='WNT_gene_list')

an_params.add_gene_list(markers = ['EGR1','RSPO2','FGFR4','WNT2','CDO1','FGF7','BMP5','BMP4',
								   'LIFR','TAGLN','ACTA2','MYL9','MYH11','ACTG2','HHIP','FOXF1',
								   'WNT5A','NGFR','SOX11','PDGFRA','PTCH1','HOXB5','HOXB4'],
					    label='explant_mesenchyme_list')

an_params.add_gene_list(markers = ['RSPO2','WNT2','CDO1','BMP5','LIFR','FGFR4','EGR1','TAGLN',
								   'ACTA2','MYL9','MYH11','ACTG2','HHIP','FOXF1','WNT5A','ELN',
								   'MGP','BGN','AGTR2','NGFR','FGF7','BMP4','SOX11','PDGFRA',
								   'PTCH1','HOXB5','HOXB4'],
						label='mesenchyme_list_updated2')

an_params.add_gene_list(markers = ['SOX9','SFTPC','NPC2','TESC','CA2','ETV5','SOX2','TP63',
								   'KRT5','KRT15','IL33','FOXL1','FOXJ1','CHGA','SYN1',
								   'MUC5AC','MUC5B','FOXA3','SCGB1A1','BPIFA1','SCGB3A2',
								   'CFTR','SFTPB','HOPX','PDPN','AQP5','AGER','ABCA3'],
					 label='lung_epithelium_markers')

## Parameters used to filter the data - Mainly used to get rid of bad cells
an_params.set_qc_params(min_cells = 0,
						min_genes = 750, # Filter out cells with fewer genes to remove dead cells
						max_genes = 3000, # Filter out cells with more genes to remove most doublets
						max_counts = 15000, # Filter out cells with more UMIs to catch a few remaining doublets
						max_mito = 5,
						doublet_detection=False) # Filter out cells with high mitochondrial gene content

## Parameters used for initial clustering analysis
an_params.set_analysis_params(n_neighbors = 50, # Size of the local neighborhood used for manifold approximation
							n_pcs = 30, # Number of principle components to use in construction of neighborhood graph
							spread = 1, # In combination with min_dist determines how clumped embedded points are
							min_dist = 0.4, # Minimum distance between points on the umap graph
							resolution = 0.225,
							clustering_choice = 'louvain',
							do_bbknn = True) # High resolution attempts to increases # of clusters identified

an_params.set_plot_params(size = 3,
						   umap_obs = ['louvain','sampleName','age'],
						   rank_grouping = ['louvain'],
						   dot_grouping = ['louvain'])

## Basic pipeline for analysis - will filter data, process, cluster, etc. and output relevant figures
an_run = SCARunner.SCARunner()
# an_run.pipe_basic(an_params,figdir,load_save='adata_save.p')

## If you find some interesting clusters that you want to "zoom in" on and recluster, you can use the following code
# Will 
# New analysis parameters for the subset of parameters
an_params.set_analysis_params(n_neighbors = 30,
						n_pcs = 20,
						spread = 1,
						min_dist = 0.4,
						resolution = 0.55,
						clustering_choice = 'louvain',
						dpt = ['louvain',['0']],
						do_bbknn = False,
						draw_force_atlas = True,
						phate=True,
						umap_init_pos='paga')

# an_run.pipe_ext(an_params, figdir=figdir,load_save='adata_save.p',label='3_extracted',extracted=['3'])

# an_run.pipe_ext(an_params, figdir=figdir,load_save='adata_save.p',label='0_1_2_4_extracted',extracted=['0','1','2','4'])

# an_run.clusters2_compare = ['all',['0','1']]
# an_run.pipe_ext(an_params, figdir=figdir,load_save='adata_save.p',label='0_1_extracted',extracted=['0','1'])

# ##### Slightly adjusting qc parameters for cleaner dpt and force atlas graphs #####
# import pickle
# run_save = pickle.load(open(''.join([figdir,'/extracted/0_1_extracted/extracted_adata_save.p']),"rb"))
# adata_postQC = run_save.adata_postQC[run_save.adata_postQC.obs['n_counts'] < 8000].copy()
# # an_run.pipe_basic(an_params, adata_filtered = adata_postQC, figdir=''.join([figdir,'/extracted/0_1_extracted_reqc_bbknn/']))
# an_params.annotation_dict = run_save.annotation_dict.copy()
# # print(run_save.annotation_dict.copy())
# for sample in an_params.sample_list:
# 	sample_name = dict([i.split(':') for i in an_params.annotation_dict[sample][1:]])['sampleName']
# 	adata_postQC_sub_sample = adata_postQC[adata_postQC.obs['sampleName']==sample_name].copy()
# 	an_run.pipe_basic(an_params, adata_filtered = adata_postQC_sub_sample, figdir=''.join([figdir,'/extracted/0_1_',sample_name,'_reqc']))

# an_run.pipe_basic(an_params, adata_filtered = adata_postQC.copy(), figdir=''.join([figdir,'/extracted/0_1_extracted_phate_bbknn/']))

# an_params.analysis_params.do_bbknn = False
# an_run.pipe_basic(an_params, adata_filtered = adata_postQC.copy(), figdir=''.join([figdir,'/extracted/0_1_extracted_phate/']))

# an_run.pipe_basic(an_params, figdir=''.join([figdir,'/extracted/0_1_extracted/']), load_save='extracted_adata_save.p')

##### Trying PHATE #####
# import phate 
# X_phate = phate.PHATE(n_jobs=-2).fit_transform(adata_postQC)
# an_run.pipe_phate(an_params, figdir=figdir)

## Parameters used to filter the data - Mainly used to get rid of bad cells
an_params.set_qc_params(min_cells = 0,
						min_genes = 750, # Filter out cells with fewer genes to remove dead cells
						max_genes = 3000, # Filter out cells with more genes to remove most doublets
						max_counts = 8000, # Filter out cells with more UMIs to catch a few remaining doublets
						max_mito = 0.05,
						doublet_detection=False) # Filter out cells with high mitochondrial gene content

##### Trying CellRank
## Parameters used for initial clustering analysis
an_params.set_analysis_params(n_neighbors = 50, # Size of the local neighborhood used for manifold approximation
							n_pcs = 30, # Number of principle components to use in construction of neighborhood graph
							spread = 1, # In combination with min_dist determines how clumped embedded points are
							min_dist = 0.4, # Minimum distance between points on the umap graph
							resolution = 0.3,
							clustering_choice = 'louvain',
							do_bbknn = True) # High resolution attempts to increases # of clusters identified

# an_run.pipe_basic(an_params,figdir=figdir,load_file_type='.loom', load_save='adata_save.p')

an_params.set_analysis_params(n_neighbors = 30,
						n_pcs = 20,
						spread = 1,
						min_dist = 0.4,
						resolution = 0.55,
						clustering_choice = 'louvain',
						# dpt = ['louvain',['0']],
						do_bbknn = False,
						draw_force_atlas = False,
						phate=False)
						# umap_init_pos='paga')

# an_run.pipe_ext(an_params, figdir=figdir,load_save='adata_save.p',label='0_1_extracted',extracted=['0','1'])

figdir = './figures_070521_cellrank/'
import pickle
run_save = pickle.load(open(''.join([figdir,'/extracted/0_1_extracted/extracted_adata_save.p']),"rb"))
adata_postQC = run_save.adata_postQC.copy()
an_params.annotation_dict = run_save.annotation_dict.copy()
# # del run_save
# # an_run.pipe_cell_rank(an_params, adata_filtered=adata_postQC, figdir=figdir)
# an_run.pipe_basic(an_params, adata_filtered=adata_postQC.copy(), figdir=''.join([figdir,'/extracted/0_2_extracted_reqc/']))

for sample in an_params.sample_list:
	sample_name = dict([i.split(':') for i in an_params.annotation_dict[sample][1:]])['sampleName']
	adata_postQC_sub_sample = adata_postQC[adata_postQC.obs['sampleName']==sample_name].copy()
	an_run.pipe_basic(an_params, adata_filtered=adata_postQC_sub_sample, figdir=''.join([figdir,'/extracted/0_1_',sample_name,'_cell_rank/']))
	# an_run.pipe_cell_rank(an_params, adata_filtered=adata_postQC_sub_sample, figdir=''.join([figdir,'/extracted/0_1_',sample_name,'_cell_rank/']))