'''
Distal Fetal Lung CELL ANALYSIS SCRIPT
By Josh Wu
13 April, 2020


Contains of Distal Lung for Renee Conway
'''
import sys
sys.path.insert(0,'C:/Users/Josh/Desktop/sca/tools')
#from scanpy_spence import *
from sca_run import *

figdir = './figures_052120/'
an_run = sca_run()

#############################################################################
## Change this to point toward your mount location for our MiStorage share ##
#############################################################################
an_run.storage_mount_point = 'Z:/'

## IDs of samples as represented in the metadata table
an_run.sample_list = ['2268-2','2598-33','2321-2','2288-4','2444-1'] # Distal Fetal Lung data

## List of interesting genes
an_run.add_gene_list(markers=['WNT2','WNT2B','WNT7B','RSPO1','RSPO3','RSPO4'],label='basic_list')
	#['EPCAM','VIM','TCF21','CDH5','CD74','NEUROD1','TAGLN','RSPO2','FGFR4','WNT2'])

an_run.add_gene_list(markers=['GPC3','CHST2','PIEZO2','ABCA8','SLC40A1','LGALS1','SLIT2','ITM2C','SLC1A5', 
							  'BNIP3L','MMP23B','PLXDC2','APOE','GJA5','SCN7A','TMEM108','S1PR1','ABCA6', 
							  'LIFR','ELN','MYLK','HHIP','MYLK','PLP1','BGN','SPARC','CD82','TYRP1', 
							  'SERPINE2','WNT5A','KCNMB1','MRGPRF','PDGFC','KCNK17','MATN2','SVIL','AGTR2', 
							  'WNT3A','PDGFRA','RSPO2','TAGLN','FGFR4','VIM','EPCAM','CD74','CDH5','WNT7B'],
							  label='extended_list')

an_run.add_gene_list(markers = ['RSPO2','FGFR4','WNT2','TAGLN','ACTA2','HHIP','PDGFRB',
								'PDGFRA','RSPO1','RSPO3','RSPO4','LGR4','LGR5','LGR6',
								'AXIN2','SOX9','SFTPC'],
					 label = 'general_list')

an_run.add_gene_list(markers = ['RSPO2','FGFR4','WNT2','TAGLN','ACTA2','HHIP','PDGFRB',
								'PDGFRA','RSPO1','RSPO3','RSPO4','LGR4','LGR6'],
					 label = 'mesenchyme_list')

an_run.add_gene_list(markers = ['SOX11','NGFR','PTCH1','HOXB5','HOXB4','BMP4','TGFB1','TGFB2',
								'TGFB3','FOXF1','EGR1','TCF21','FN1','CDO1','BOC','PTCH1',
								'SELL','FGF7','BMP5','MYL9','MYH11','ACTG2','PDGFRB','EGFR',
								'NOTCH2','BMPR1A','BMPR2'],
					 label = 'extracted_list')

an_run.add_gene_list(markers = ['EPCAM','CDH1','KRT18','KRT8','CLDN6','VIM','POSTN','DCN','TCF21',
								 'COL1A2','COL3A1','COL1A2','TAGLN','ACTA2','PDGFRB','CD151','PTPRC','CD37',
								 'CORO1A','LCP1','CD53','LAPTM5','CDH5','CLDN5','ESAM','KDR','FLT1','S100B',
								 'ELAVL4','TUBB2B','STMN2','ASCL1','NNAT','GRP','MPZ'],
					  label= 'gene_list_2')

an_run.add_gene_list(markers = ['FGF10','FOS','MKI67','TOP2A','CCND1','TYMS','DUT','CDK1','H2AFZ','STMN1','FAS','CASP3'],
					  label='mesenchyme_extracted_list')

## Parameters used to filter the data - Mainly used to get rid of bad cells
an_run.set_filter_params(min_cells = 0,
						min_genes = 750, # Filter out cells with fewer genes to remove dead cells
						max_genes = 3000, # Filter out cells with more genes to remove most doublets
						max_counts = 15000, # Filter out cells with more UMIs to catch a few remaining doublets
						max_mito = 0.05) # Filter out cells with high mitochondrial gene content

## Parameters used for initial clustering analysis
an_run.set_analysis_params(n_neighbors = 30, # Size of the local neighborhood used for manifold approximation
							n_pcs = 20, # Number of principle components to use in construction of neighborhood graph
							spread = 1, # In combination with min_dist determines how clumped embedded points are
							min_dist = 0.4, # Minimum distance between points on the umap graph
							resolution = 0.25,
							remove_batch_effects = True) # High resolution attempts to increases # of clusters identified

an_run.set_plot_params(size = 3,
					   umap_obs = ['louvain','sampleName','age'],
					   exp_grouping = ['louvain'])

## Basic pipeline for analysis - will filter data, process, cluster, etc. and output relevant figures
# an_run.pipe_basic(figdir,load_save='adata_save.p')

## If you find some interesting clusters that you want to "zoom in" on and recluster, you can use the following code
# Will 
# New analysis parameters for the subset of parameters
analysis_params_ext = dict(n_neighbors = 15,
						n_pcs = 10,
						spread = 1,
						min_dist = 0.4,
						resolution = 0.4,
						remove_batch_effects = True)

# an_run.pipe_ext(analysis_params_ext, figdir=figdir,load_save='adata_save.p',label='0_1_extracted',extracted=['0','1'])

an_run.pipe_ext(analysis_params_ext, figdir=figdir,load_save='adata_save.p',label='0_1_2_4_extracted',extracted=['0','1','2','4'])

an_run.clusters2_compare = ['all',['0','1','2']]
an_run.pipe_ext(analysis_params_ext, figdir=figdir,load_save='adata_save.p',label='0_1_2_extracted',extracted=['0','1','2'])#, load_save='adata_save.p')

