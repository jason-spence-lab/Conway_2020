'''
BASIC SINGLE CELL ANALYSIS SCRIPT
by Josh Wu
4 June, 2019

Relies heavily on the Scanpy Python module developed by the Theis Lab
Read more about Scanpy at https://scanpy.readthedocs.io/en/latest/index.html

Contains analysis of cultured lung samples for Renee Conway

In progress ---
Moving to encapsulate parameters and relevant functions using class sca_set()
'''
import sys
sys.path.insert(0,'C:/Users/Josh/Desktop/sca_run/sca_pipes')
#from scanpy_spence import *
from classes.sca_params import *
import SCARunner as SCARunner

#from tools.pipelines import *
def run_analysis():
	figdir = './figures_071721/'
	an_params = sca_params()
	#############################################################################
	## Change this to point toward your mount location for our MiStorage share ##
	#############################################################################
	an_params.storage_mount_point = 'Z:/'

	# ## List of interesting genes

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

	an_params.add_gene_list(markers = ['SOX9','SFTPC','NPC2','TESC','CA2','ETV5','SOX2','TP63',
									   'KRT5','KRT15','IL33','FOXL1','FOXJ1','CHGA','SYN1',
									   'MUC5AC','MUC5B','FOXA3','SCGB1A1','BPIFA1','SCGB3A2',
									   'CFTR','SFTPB','HOPX','PDPN','AQP5','AGER','ABCA3'],
						 label='lung_epithelium_markers')

	# cell_score_lists = ['BTP_short','BTP_long','basal_score',
 #   					    'budtipadj_score','club_score','goblet_score',
 #   					    'int_ciliated_score','multiciliated_score',
 #   					    'neuroendocrine_score','multiciliated_pre_score',
 #   					    'secretory_progenitor_score','submucosa_score',
 #   					    'submucosa_basal_score','basal_full','budtipadj_full',
 #   					    'club_full','goblet_full','int_ciliated_full','multiciliated_full',
	# 					'neuroendocrine_full','multiciliated_pre_full',
	# 					'secretory_progenitor_full','submucosa_full',
	# 					'submucosa_basal_full','WNT_CellScore']

	cell_score_lists = ['Multiciliated_CellScore_LogFold','SecretoryProgenitor_CellScore_LogFold',
						'SubmucosaBasal_CellScore_LogFold','MulticiliatedPrecursor_CellScore_LogFold',
						'BudTipProgenitor_CellScore_LogFold','BudTipAdjacent_CellScore_LogFold',
						'IntermediateCiliated_CellScore_LogFold','BasalCell_CellScore_LogFold',
						'Submucosa_CellScore_LogFold','GobletCell_CellScore_LogFold','ClubCell_CellScore_LogFold']
	
	for score_list in cell_score_lists:
   		an_params.add_gene_list(load_file=score_list, label=score_list, cell_score_list='Only')

	## Parameters used for initial clustering analysis
	an_params.set_analysis_params(n_neighbors = 30, # Size of the local neighborhood used for manifold approximation
							   n_pcs = 20, # Number of principle components to use in construction of neighborhood graph
							   spread = 1, # In combination with min_dist determines how clumped embedded points are
							   min_dist = 0.4, # Minimum distance between points on the umap graph.
							   clustering_choice = 'louvain',
							   resolution = 0.4) # High resolution attempts to increases # of clusters identified

	## Parameters used to filter the data - Mainly used to get rid of bad cells
	an_params.set_qc_params(min_cells = 0, # Filter out cells 
							 min_genes = 750, # Filter out cells with fewer genes to remove dead cells
							 max_genes = 10000, # Filter out cells with more genes to remove most doublets
							 max_counts = 50000, # Filter out cells with more UMIs to catch a few remaining doublets
							 max_mito = 0.1) # Filter out cells with high mitochondrial gene content

	an_params.set_plot_params(umap_obs = ['louvain','sampleName'],
						   dot_grouping = ['louvain'],
						   rank_grouping = ['louvain','sampleName'],
						   size=5)#,
						   # final_quality=True)
	an_run = SCARunner.SCARunner()

	an_params.sample_list=['283-1','283-2','283-3']
	# an_run.pipe_basic(an_params,''.join([figdir,'all/']))#,load_save='adata_save.p')

	an_params.sample_list=['283-2','283-3']
	# an_run.pipe_basic(an_params,''.join([figdir,'control_+_lgd_ecd/']))#,load_save='adata_save.p')

	# an_params.plot_params.rank_grouping = ['louvain']

	# an_params.sample_list=['283-1']
	# an_run.pipe_basic(an_params,''.join([figdir,'no_virus/']))#,load_save='adata_save.p')

	# an_params.sample_list=['283-2']
	# an_run.pipe_basic(an_params,''.join([figdir,'control_virus/']))#,load_save='adata_save.p')

	# an_params.sample_list=['283-3']
	# an_run.pipe_basic(an_params,''.join([figdir,'lgd_ecd/']))#,load_save='adata_save.p')
	
	## If you find some interesting clusters that you want to "zoom in" on and recluster, you can use the following code
	# New analysis parameters for the subset of parameters
	an_params.set_analysis_params(n_neighbors = 20,
								n_pcs = 15,
								spread = 1,
								min_dist = 0.4,
								resolution = 0.4,
								clustering_choice = 'louvain',)

	an_params.plot_params.size=20
	an_params.plot_params.rank_grouping = ['louvain','sampleName']
	# an_run.pipe_ext(an_params, figdir=''.join([figdir,'all/']), label='epithelium_1_8',
	#  			extracted=['1','8'], load_save='adata_save.p')

	# an_run.pipe_ext(an_params, figdir=''.join([figdir,'control_+_lgd_ecd/']), label='epithelium_1_8',
	#  			extracted=['1','8'], load_save='adata_save.p')

	an_run.pipe_ext(an_params, figdir=''.join([figdir, 'control_+_lgd_ecd/extracted/epithelium_1_8/']),
					extracted=['0'], label='c0', load_save='extracted_adata_save.p')

	# ## Violin plots for filtering parameters pre and post
	# sc.pl.violin(an_params.adata, an_params.gene_dict['epithelial_markers']['markers'],groupby='louvain',
	# 			 jitter=0.4, multi_panel=True, save='_epithelial_louvain.png', show=False)

	# sc.pl.violin(an_params.adata, an_params.gene_dict['epithelial_markers']['markers'],groupby='sampleName',
	# 			 jitter=0.4, multi_panel=True, save='_epithelial_sample.png', show=False)

	# an_params.plot_params.rank_grouping = ['louvain']
	# an_run.pipe_ext(an_params, figdir=''.join([figdir,'no_virus/']), label='epithelium_3_5_6',
	#  			extracted=['3','5','6'], load_save='adata_save.p')

	# an_run.pipe_ext(an_params, figdir=''.join([figdir,'control_virus/']), label='epithelium_2',
	#  			extracted=['2'], load_save='adata_save.p')

	# an_run.pipe_ext(an_params, figdir=''.join([figdir,'lgd_ecd/']), label='epithelium_2_7',
	#  			extracted=['2','7'], load_save='adata_save.p')

							
if __name__ == "__main__":
	run_analysis()