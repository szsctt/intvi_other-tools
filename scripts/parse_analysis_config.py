import pdb
import itertools
import pandas as pd
from os import path
from glob import glob

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
			analysis_conditions += make_polyidus_rows(config, dataset)
			
		if 'vifi_params' in config[dataset]:
			analysis_conditions += make_vifi_rows(config, dataset)
			
		if 'verse_params' in config[dataset]:
			analysis_conditions += make_verse_rows(config, dataset)

	
	# make data frame 
	return pd.DataFrame(analysis_conditions, columns = column_names)
			
def make_verse_rows(config, dataset):
	#### parameters for verse ####
	rows = []
	
	# paramters should be lists, so that we can do combinations of all
	assert hasattr(config[dataset]['verse_params']['detection_mode'], '__iter__')
	assert hasattr(config[dataset]['verse_params']['flank_region_size'], '__iter__')
	assert hasattr(config[dataset]['verse_params']['sensitivity_level'], '__iter__')
	assert hasattr(config[dataset]['verse_params']['min_contig_length'], '__iter__')
	assert hasattr(config[dataset]['verse_params']['blastn_evalue_thrd'], '__iter__')	
	assert hasattr(config[dataset]['verse_params']['similarity_thrd'], '__iter__')
	assert hasattr(config[dataset]['verse_params']['chop_read_length'], '__iter__')
	assert hasattr(config[dataset]['verse_params']['minIdentity'], '__iter__')

 
	# each combination of these are a unique 'analysis condition' for our pipeline
	i = 0 
	for host, virus, detection_mode, flank_region_size, sensitivity_level, min_contig_length, blastn_evalue_thrd, similarity_thrd, chop_read_length, minIdentity in itertools.product(
																		config[dataset]['analysis_hosts'].keys(),
																		config[dataset]['analysis_viruses'].keys(),
																		config[dataset]['verse_params']['detection_mode'], 
																		config[dataset]['verse_params']['flank_region_size'], 
																		config[dataset]['verse_params']['sensitivity_level'],
																		config[dataset]['verse_params']['min_contig_length'],
																		config[dataset]['verse_params']['blastn_evalue_thrd'],
																		config[dataset]['verse_params']['similarity_thrd'],
																		config[dataset]['verse_params']['chop_read_length'],
																		config[dataset]['verse_params']['minIdentity'],
																		):
		
		condition = f"{dataset}_verse{i}"
		rows.append({
				'experiment' : dataset,
				'host' 			: host,
				'host_fasta': config[dataset]['analysis_hosts'][host],
				'virus'     : virus,
				'virus_fasta': config[dataset]['analysis_viruses'][virus],
				'analysis_condition': condition,
				'detection_mode' : detection_mode,
				'flank_region_size' : flank_region_size,
				'sensitivity_level' : sensitivity_level,
				'min_contig_length' : min_contig_length,
				'blastn_evalue_thrd': blastn_evalue_thrd,
				'similarity_thrd'   : similarity_thrd,
				'chop_read_length'  : chop_read_length,
				'minIdentity'       : minIdentity,
				'tool'			 : 'verse',	
			})
		i += 1
	return add_read_info(config, dataset, rows)
	

def make_vifi_rows(config, dataset):
	#### parameters for vifi ####
	i = 0
	rows = []
	# make sure the required information about the host genome has been provided
	host_file_keys = ('mappability', 'mappability_exclude', 'genes', 'exons', 'oncogenes', 'centromeres',
												'conserved_regions', 'segdup')
	hosts_to_use = []
	for host in config[dataset]['analysis_hosts'].keys():
		# check that this host is in host_info
		if host not in config[dataset]['vifi_params']['host_info']:
			print(f"host_info not provided: skipping ViFi for {host}")
			continue
		# check that all necessary files have been specified
		if not all([key in config[dataset]['vifi_params']['host_info'][host] for key in host_file_keys]):
			print(f"one or more of the required files ({host_file_keys}) for host {host} is not specfied: skipping ViFi for host {host}")
			continue
		hosts_to_use.append(host)
									

	for host, virus in itertools.product(hosts_to_use, config[dataset]['analysis_viruses'].keys()):
		condition = f"{dataset}_vifi{i}"
		rows.append({
				'experiment' : dataset,
				'host' 			 : host,
				'host_fasta' : config[dataset]['analysis_hosts'][host],
				'host_mappability' : config[dataset]['vifi_params']['host_info'][host]['mappability'],
				'host_mappability_exclude' : config[dataset]['vifi_params']['host_info'][host]['mappability_exclude'],				
				'host_genes' : config[dataset]['vifi_params']['host_info'][host]['genes'],				
				'host_exons' : config[dataset]['vifi_params']['host_info'][host]['exons'],				
				'host_oncogenes' : config[dataset]['vifi_params']['host_info'][host]['oncogenes'],	
				'host_centromeres' : config[dataset]['vifi_params']['host_info'][host]['centromeres'],	
				'host_conserved_regions' : config[dataset]['vifi_params']['host_info'][host]['conserved_regions'],		
				'host_segdup' : config[dataset]['vifi_params']['host_info'][host]['segdup'],																
				'virus'      : virus,		
				'virus_fasta': config[dataset]['analysis_viruses'][virus],	
				'analysis_condition': condition,
				'tool'			 : 'vifi',			
				})
		i += 1
	return add_read_info(config, dataset, rows)


def make_polyidus_rows(config, dataset):
	#### parameters for polyidus ####
	rows = []
	i = 0
	
	# are we trying multiple aligners?
	if 'aligner' in config[dataset]['polyidus_params']:
		assert hasattr(config[dataset]['polyidus_params']['aligner'], '__iter__')
		aligners = config[dataset]['polyidus_params']['aligner']
	else:
		aligners = ['bowtie2']
				
	for host, virus, aligner in itertools.product(
																		config[dataset]['analysis_hosts'].keys(),
																		config[dataset]['analysis_viruses'].keys(),
																		aligners):
		# give this analysis condition a name
		condition = f"{dataset}_polyidus{i}"

		rows.append({
					'experiment' : dataset,
					'host' 			 : host,
					'host_fasta' : config[dataset]['analysis_hosts'][host],
					'virus'      : virus,
					'virus_fasta': config[dataset]['analysis_viruses'][virus],
					'analysis_condition': condition,
					'aligner'		 : aligner,
					'merge'			 : 0,
					'tool'			 : 'polyidus',
					})	
		i += 1
	return add_read_info(config, dataset, rows)

def get_bool_value_from_config(config, dataset, key, default):
	if key not in config[dataset]:
		return default
		
	if config[dataset][key] is True:
		return 1
	elif config[dataset][key] is False:
		return 0
	else:
		raise ValueError(f"Boolean value for {key} in dataset {dataset} is neither True nor False.  Please specify one or the other")
		
def add_read_info(config, dataset, rows):
		# add 'read_folder', 'R1_suffix', 'R2_suffix', 'outdir', 'bam_suffix', if they are defined
 		read_folder = config[dataset].get('read_directory')
 		outdir = config[dataset].get('read_directory')
 		R1_suffix = config[dataset].get('R1_suffix')
 		R2_suffix = config[dataset].get('R2_suffix')
 		bam_suffix = config[dataset].get('bam_suffix')
 		
 		for row in rows:
 			row['read_folder'] = read_folder
 			row['outdir'] = outdir
 			row['R1_suffix'] = R1_suffix
 			row['R2_suffix'] = R2_suffix
 			row['bam_suffix'] = bam_suffix
 			
 		return rows
	
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
		
		
		
		
