### snakemake workflow to analyse paired-end reads using various viral integration tools

import pdb
import itertools
import pandas as pd

from scripts.parse_analysis_config_single import parse_analysis_config, get_samples
from scripts.input_functions import *

verse_threads = 8

#####################################################
############## analysis parameters ##################
#####################################################

# parse config file to get desired analyses
analysis_df = parse_analysis_config(config)

# get samples for each dataset
sample_dict = get_samples(config)

# make dict with host/virus names as keys and paths to fasta files as values
ref_names = {name: path for name, path in zip(list(analysis_df.host) + list(analysis_df.virus), 	
				list(analysis_df.host_fasta) + list(analysis_df.virus_fasta))} 
# check that each host/virus name is unique
if len(ref_names.keys()) != len(set(analysis_df.host)) + len(set(analysis_df.virus)):
	raise ValueError("Each host and virus reference name must be unique")				

#####################################################
################## wildcards ########################
#####################################################

all_samples = set()
for samples in sample_dict.values():
	all_samples |= set(samples)

wildcard_constraints:
	outpath = "|".join(set(analysis_df['outdir'])),
	dset = "|".join(set(analysis_df['experiment'])),
	host = "|".join(set(analysis_df['host'])),
	virus = "|".join(set(analysis_df['virus'])),	
	sample = "|".join(set(all_samples))


#####################################################
################### target files ####################
#####################################################

analysis_summaries = expand("{outdir}/{experiment}/analysis_conditions.tsv", 
			zip,
			exp = analysis_df.experiment, 
			outpath = analysis_df.outdir
			)


# targets for other tools
other_tool_targets = set()
for i, row in analysis_df.iterrows():

	samples = sample_dict[row['experiment']]
	
	if row['tool'] == 'polyidus':
		other_tool_targets |= set(expand("{outpath}/{dset}/polyidus/{host}.{virus}.{samp}/results/exactHpvIntegrations.tsv", 
																	outpath = row['outdir'],
																	dset = row['experiment'],
																	analysis_condition = row['analysis_condition'],
																	host = row['host'],
																	virus = row['virus'],
																	samp = samples
														))
	elif row['tool'] == 'seeksv':
		other_tool_targets |= set(expand("{outpath}/{dset}/seeksv/ints/{samp}.{host}.{virus}.integrations.txt", 
																	outpath = row['outdir'],
																	dset = row['experiment'],
																	analysis_condition = row['analysis_condition'],
																	host = row['host'],
																	virus = row['virus'],
																	samp = samples
														))
	elif row['tool'] == 'verse':
		other_tool_targets |= set(expand("{outpath}/{dset}/verse/{samp}.{host}.{virus}/integration-sites.txt", 
																	outpath = row['outdir'],
																	dset = row['experiment'],
																	analysis_condition = row['analysis_condition'],
																	host = row['host'],
																	virus = row['virus'],
																	samp = samples
														))				
	elif row['tool'] == 'vifi':
		other_tool_targets |= set(expand("{outpath}/{dset}/vifi/{samp}.{host}.{virus}/output.clusters.txt", 
																	outpath = row['outdir'],
																	dset = row['experiment'],
																	analysis_condition = row['analysis_condition'],
																	host = row['host'],
																	virus = row['virus'],
																	samp = samples
														))			
	elif row['tool'] == 'vseq_toolkit':
		other_tool_targets |= set(expand("{outpath}/{dset}/vseq_toolkit/{samp}.{host}.{virus}/ISGenomeVector.csv", 
																	outpath = row['outdir'],
																	dset = row['experiment'],
																	analysis_condition = row['analysis_condition'],
																	host = row['host'],
																	virus = row['virus'],
																	samp = samples
														))																		


rule all: 
	input:
		set(analysis_summaries),
		other_tool_targets
	
include: "snakemake_rules/trim.smk"

include: "snakemake_rules/polyidus.smk"

include: "snakemake_rules/seeksv.smk"

include: "snakemake_rules/verse.smk"

include: "snakemake_rules/vifi.smk"

include: "snakemake_rules/vseq.smk"
