import pdb
import itertools
import pandas as pd
from os import path
from glob import glob

# defaults
polyidus_default_aligner = 'bowtie2'
polyidus_default_trim = False

vifi_default_trim = False

seeksv_default_trim = False
seeksv_default_dedup = False

verse_default_trim = False
verse_default_detection_mode = 'sensitive'
verse_default_flank_region_size = 4000
verse_default_sensitivity_level = 1
verse_default_min_contig_length = 300
verse_default_blastn_evalue_thrd = 0.05
verse_default_similarity_thrd = 0.8
verse_default_chop_read_length = 25
verse_default_minIdentity = 80


def parse_analysis_config(config):
	
	# global parameters (defaults for this dataset) may be specified using a 'global' key
	# at the top of the config file.  These will be applied to any missing keys in the remaining 
	# datasets in the config file.  If a key is specified in both 'global' and a dataset,
	# the value specified in the dataset will be used
	
	if 'global' in config:		
		# get default (global) options
		default = config.pop('global')
		for dataset in config:
			for key in default:
				if key not in config[dataset]:
					config[dataset][key] = default[key]
					
	# get unique analysis conditions - these are combinations of the analysis parameters that
	# can be set in our pipeline (merging, de-duplicaiton, bwa mem prarameters, etc), or
	# the tool to be used in analysis

	column_names = ('experiment', 'exp', 'analysis_condition', 
									'tool', 'host', 'host_fasta',
									'virus', 'virus_fasta', 'read_folder', 
									'R1_suffix', 'R2_suffix', 'outdir', 'bam_suffix',
									'adapter_1', 'adapter_2', 'merge', 'trim', 'dedup',
									'host_mappability', 'host_mappability_exclude', 
									'host_genes', 'host_exons', 'host_oncogenes', 
									'host_centromeres', 'host_conserved_regions',
									'host_segdup', 'detection_mode', 
									'flank_region_size', 'sensitivity_level', 
									'min_contig_length', 'blastn_evalue_thrd', 
									'similarity_thrd', 
									'chop_read_length', 'minIdentity')		

	analysis_conditions = []
	
	for dataset in config.keys():
		
		if 'polyidus_params' in config[dataset]:
			analysis_conditions.append(make_polyidus_row(config, dataset))
			
		if 'vifi_params' in config[dataset]:
			analysis_conditions.append(make_vifi_row(config, dataset))
			
		if 'verse_params' in config[dataset]:
			analysis_conditions.append(make_verse_row(config, dataset))

		if 'seeksv_params' in config[dataset]:
			analysis_conditions.append(make_seeksv_row(config, dataset))		
	
	# make data frame 
	return pd.DataFrame(analysis_conditions, columns = column_names)
			
def make_verse_row(config, dataset):
	#### parameters for verse ####

	trim = get_bool(config[dataset]['verse_params'], 'trim', verse_default_trim, 'verse_params')	
	detection_mode = get_entry_with_default(config[dataset]['verse_params'], 
																							'detection_mode', verse_default_detection_mode, 'verse_params')	
	flank_region_size = get_entry_with_default(config[dataset]['verse_params'], 
																							'flank_region_size', 
																							verse_default_flank_region_size , 'verse_params')	
	sensitivity_level = get_entry_with_default(config[dataset]['verse_params'], 
																							'sensitivity_level', 
																							verse_default_sensitivity_level, 'verse_params')		
	min_contig_length = get_entry_with_default(config[dataset]['verse_params'], 
																							'min_contig_length', 
																							verse_default_min_contig_length, 'verse_params')
	blastn_evalue_thrd = get_entry_with_default(config[dataset]['verse_params'], 
																							'blastn_evalue_thrd', 
																							verse_default_blastn_evalue_thrd, 'verse_params')
	similarity_thrd = get_entry_with_default(config[dataset]['verse_params'], 
																							'similarity_thrd', 
																							verse_default_similarity_thrd, 'verse_params')
	chop_read_length = get_entry_with_default(config[dataset]['verse_params'], 
																							'chop_read_length', 
																							verse_default_chop_read_length, 'verse_params')
	minIdentity = get_entry_with_default(config[dataset]['verse_params'], 
																							'minIdentity', 
																							verse_default_minIdentity, 'verse_params')

 
	row = {
				'detection_mode' : detection_mode,
				'flank_region_size' : flank_region_size,
				'sensitivity_level' : sensitivity_level,
				'min_contig_length' : min_contig_length,
				'blastn_evalue_thrd': blastn_evalue_thrd,
				'similarity_thrd'   : similarity_thrd,
				'chop_read_length'  : chop_read_length,
				'minIdentity'       : minIdentity,
				'trim'              : trim,
				'merge'             : 0,
				'dedup'             : 0,
				'tool'			 : 'verse',	
			}
	return add_common_info(config, dataset, row)
	

def make_vifi_row(config, dataset):
	#### parameters for vifi ####
	# make sure the required information about the host genome has been provided
	host_file_keys = ('mappability', 'mappability_exclude', 'genes', 'exons', 'oncogenes', 'centromeres',
												'conserved_regions', 'segdup')
	
	assert len(config[dataset]['analysis_host'].keys()) == 1
	host = list(config[dataset]['analysis_host'].keys())[0]
	assert host in config[dataset]['vifi_params']['host_info']
	
	# check that all necessary files have been specified
	if not all([key in config[dataset]['vifi_params']['host_info'][host] for key in host_file_keys]):
		print(f"one or more of the required files ({host_file_keys}) for host {host} is not specfied")
		return {}
		
	trim = get_bool(config[dataset]['polyidus_params'], 'trim', polyidus_default_trim, 'polyidus_params')						

	row = {
				'host_mappability' : config[dataset]['vifi_params']['host_info'][host]['mappability'],
				'host_mappability_exclude' : config[dataset]['vifi_params']['host_info'][host]['mappability_exclude'],				
				'host_genes' : config[dataset]['vifi_params']['host_info'][host]['genes'],				
				'host_exons' : config[dataset]['vifi_params']['host_info'][host]['exons'],				
				'host_oncogenes' : config[dataset]['vifi_params']['host_info'][host]['oncogenes'],	
				'host_centromeres' : config[dataset]['vifi_params']['host_info'][host]['centromeres'],	
				'host_conserved_regions' : config[dataset]['vifi_params']['host_info'][host]['conserved_regions'],		
				'host_segdup' : config[dataset]['vifi_params']['host_info'][host]['segdup'],																
				'tool'			 : 'vifi',	
				'trim'			 : trim,
				'merge'			 : 0,
				'dedup'      : 0,	
				}

	return add_common_info(config, dataset, row)


def make_polyidus_row(config, dataset):
	#### parameters for polyidus ####
	rows = []
	i = 0

	trim = get_bool(config[dataset]['polyidus_params'], 'trim', polyidus_default_trim, 'polyidus_params')	
	aligner = get_entry_with_default(config[dataset]['polyidus_params'], 'aligner', polyidus_default_aligner, 'polyidus_params')

	row = {
				'aligner'		 : aligner,
				'trim'       : trim,
				'dedup'      : 0,
				'tool'			 : 'polyidus',
				}

	return add_common_info(config, dataset, row)
	
def make_seeksv_row(config, dataset):
	#### parameters for polyidus ####
	
	trim = get_bool(config[dataset]['seeksv_params'], 'trim', seeksv_default_trim, 'seeksv_params')
	dedup = get_bool(config[dataset]['seeksv_params'], 'dedup', seeksv_default_dedup, 'seeksv_params')
		
	row  = {
				'trim'			 : trim,
				'dedup'      : dedup,
				'tool'			 : 'seeksv',
					}

	return add_common_info(config, dataset, row)


def get_bool_value_from_config(config, dataset, key, default):
	if key not in config[dataset]:
		return default
		
	if config[dataset][key] is True:
		return 1
	elif config[dataset][key] is False:
		return 0
	else:
		raise ValueError(f"Boolean value for {key} in dataset {dataset} is neither True nor False.  Please specify one or the other")
		
def add_common_info(config, dataset, row):
	# add 'read_folder', 'R1_suffix', 'R2_suffix', 'outdir', 'bam_suffix', if they are defined

	row['read_folder'] = path.normpath(config[dataset].get('read_directory'))
	row['outdir'] = path.normpath(config[dataset].get('out_directory'))
	row['R1_suffix'] = config[dataset].get('R1_suffix')
	row['R2_suffix'] = config[dataset].get('R2_suffix')
	row['bam_suffix'] = config[dataset].get('bam_suffix')
	row['adapter_1'] = config[dataset].get('read1-adapt')
	row['adapter_2'] = config[dataset].get('read2-adapt')
	row['experiment'] = dataset

	# make sure there's only one host and virus
	assert len(config[dataset]['analysis_host'].keys()) == 1
	assert len(config[dataset]['analysis_virus'].keys()) == 1	
	
	# add host and virus info
	row['host'] = list(config[dataset]['analysis_host'].keys())[0]
	row['virus'] =list(config[dataset]['analysis_virus'].keys())[0]
		
	row['host_fasta'] = config[dataset]['analysis_host'][row['host']]
	row['virus_fasta'] = config[dataset]['analysis_virus'][row['virus']]
 	
	return row
	
def get_samples(config):
	samples = {}
	## return a dict with datasets as key and sample names as values
	for dataset in config:
		read_dir = path.normpath(config[dataset]['read_directory'])
		samps = glob(path.join(read_dir, f"*{config[dataset]['R1_suffix']}"))
		# strip directory and suffix
		samps = [path.basename(samp)[:-len(config[dataset]['R1_suffix'])] for samp in samps]
		
		samples[dataset] = samps
		
	return samples
		
def get_bool(parent_dict, key, default, name):
	param = get_entry_with_default(parent_dict, key, default, name)
	assert isinstance(param, bool)
	return 1 if param is True else 0
	
def get_entry_with_default(parent_dict, key, default, name):
	if key not in parent_dict:
		print(f"paramter {key} not specified for {name}: using default {default}")
		return default
		
	return parent_dict[key]

