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
sys.path.insert(0,'C:/Users/Josh/Desktop/sca_run/sca_run')

from sca_run import *
import csv
#from tools.pipelines import *
def run_analysis():
	figdir = './figures_093020/'
	an_run = sca_run()
	#############################################################################
	## Change this to point toward your mount location for our MiStorage share ##
	#############################################################################
	an_run.storage_mount_point = 'Z:/'

	# ## List of interesting genes

	an_run.add_gene_list(markers = ['EPCAM','KRT18','KRT8','CLDN6','VIM','POSTN','DCN',
									'TCF21','COL1A2','TAGLN','ACTA2','PDGFRB','PTPRC',
									'CD37','CORO1A','LCP1','CD53','LAPTM5','CDH5','CLDN5',
									'ESAM','KDR','FLT1','S100B','ELAVL4','TUBB2B','STMN2',
									'NNAT','GRP','MPZ'], label='general_cell_types')

	an_run.add_gene_list(markers = ['EGR1','RSPO2','FGFR4','WNT2','CDO1','FGF7','BMP5',
									'BMP4','MYL9','MYH11','ACTG2','HHIP','FOXF1','WNT5A',
									'MKI67','TOP2A','CDK1','TYMS','H2AFZ','STMN1','SOX11',
									'NGFR','PDGFRA','PTCH1','HOXB5','HOXB4'],
						 label='mesenchyme_markers')

	an_run.add_gene_list(markers = ['RSPO1','RSPO3','RSPO4','LGR4','LGR5','LGR6','AXIN2'],
						 label='RSPO_list')

	an_run.add_gene_list(markers = ['SOX9','SFTPC','NPC2','TESC','CA2','ETV5','SOX2','TP63',
									'KRT5','KRT15','IL33','FOXI1','FOXJ1','CHGA','SYN1',
									'MUC5AC','MUC5B','FOXA3','SCGB1A1','BPIFA1',
									'SCGB3A2','CFTR','SFTPB','HOPX','PDPN','AQP5',
									'AGER','ABCA3'],
						 label='epithelial_markers')

	## Parameters used for initial clustering analysis
	an_run.set_analysis_params(n_neighbors = 20, # Size of the local neighborhood used for manifold approximation
							   n_pcs = 15, # Number of principle components to use in construction of neighborhood graph
							   spread = 1, # In combination with min_dist determines how clumped embedded points are
							   min_dist = 0.4, # Minimum distance between points on the umap graph
							   resolution = 0.4, # High resolution attempts to increases # of clusters identified
							   cell_score_lists = ['BTP_short','BTP_long','basal_score',
							   					   'budtipadj_score','club_score','goblet_score',
							   					   'int_ciliated_score','multiciliated_score',
							   					   'neuroendocrine_score','multiciliated_pre_score',
							   					   'secretory_progenitor_score','submucosa_score',
							   					   'submucosa_basal_score'])

	## Parameters used to filter the data - Mainly used to get rid of bad cells
	an_run.set_filter_params(min_cells = 0, # Filter out cells 
							 min_genes = 500, # Filter out cells with fewer genes to remove dead cells
							 max_genes = 10000, # Filter out cells with more genes to remove most doublets
							 max_counts = 50000, # Filter out cells with more UMIs to catch a few remaining doublets
							 max_mito = 0.1) # Filter out cells with high mitochondrial gene content

	an_run.set_plot_params(umap_obs = ['louvain','sampleName'],
						   exp_grouping = ['louvain'],
						   rank_grouping = ['louvain','sampleName'],
						   size=5)#,
						   # final_quality=True)

	# an_run.sample_list=['283-1','283-2','283-3']
	# an_run.pipe_basic(''.join([figdir,'all/']))#,load_save='adata_save.p')

	# an_run.sample_list=['283-2','283-3']
	# an_run.pipe_basic(''.join([figdir,'control_+_lgd_ecd/']))#,load_save='adata_save.p')

	# an_run.rank_grouping = ['louvain']

	# an_run.sample_list=['283-1']
	# an_run.pipe_basic(''.join([figdir,'no_virus/']))#,load_save='adata_save.p')

	# an_run.sample_list=['283-2']
	# an_run.pipe_basic(''.join([figdir,'control_virus/']))#,load_save='adata_save.p')

	# an_run.sample_list=['283-3']
	# an_run.pipe_basic(''.join([figdir,'lgd_ecd/']))#,load_save='adata_save.p')
	
	## If you find some interesting clusters that you want to "zoom in" on and recluster, you can use the following code
	# New analysis parameters for the subset of parameters
	analysis_params_ext = dict(n_neighbors = 15,
							n_pcs = 11,
							spread = 1,
							min_dist = 0.4,
							resolution = 0.4,
							cell_score_lists = ['BTP_short','BTP_long','basal_score',
						   					    'budtipadj_score','club_score','goblet_score',
						   					    'int_ciliated_score','multiciliated_score',
						   					    'neuroendocrine_score','multiciliated_pre_score',
						   					    'secretory_progenitor_score','submucosa_score',
						   					    'submucosa_basal_score','basal_full','budtipadj_full',
						   					    'club_full','goblet_full','int_ciliated_full','multiciliated_full',
												'neuroendocrine_full','multiciliated_pre_full',
												'secretory_progenitor_full','submucosa_full',
												'submucosa_basal_full'])

	an_run.size=20
	an_run.rank_grouping = ['louvain','sampleName']
	# an_run.pipe_ext(analysis_params_ext, figdir=''.join([figdir,'all/']), label='epithelium',
	#  			extracted=['2'], load_save='adata_save.p')

	# an_run.pipe_ext(analysis_params_ext, figdir=''.join([figdir,'control_+_lgd_ecd/']), label='epithelium',
	#  			extracted=['1'], load_save='adata_save.p')

	an_run.pipe_ext(analysis_params_ext, figdir=''.join([figdir, 'control_+_lgd_ecd/extracted/epithelium/']),
					extracted=['0'], load_save='adata_save.p')

	# ## Violin plots for filtering parameters pre and post
	sc.pl.violin(an_run.adata, an_run.gene_dict['epithelial_markers']['markers'],groupby='louvain',
				 jitter=0.4, multi_panel=True, save='_epithelial_louvain.png', show=False)

	sc.pl.violin(an_run.adata, an_run.gene_dict['epithelial_markers']['markers'],groupby='sampleName',
				 jitter=0.4, multi_panel=True, save='_epithelial_sample.png', show=False)

	# an_run.rank_grouping = ['louvain']
	# an_run.pipe_ext(analysis_params_ext, figdir=''.join([figdir,'no_virus/']), label='epithelium',
	#  			extracted=['3','5','6'], load_save='adata_save.p')

	# an_run.pipe_ext(analysis_params_ext, figdir=''.join([figdir,'control_virus/']), label='epithelium',
	#  			extracted=['2'], load_save='adata_save.p')

	# an_run.pipe_ext(analysis_params_ext, figdir=''.join([figdir,'lgd_ecd/']), label='epithelium',
	#  			extracted=['2'], load_save='adata_save.p')

							
if __name__ == "__main__":
	run_analysis()