### snakemake workflow to analyse paired-end reads using various viral integration tools

import pdb
import itertools
import pandas as pd

from scripts.parse_analysis_config import parse_analysis_config, get_samples
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
################## ruleorder ########################
#####################################################

# to resolve conflicts between polyidus and combine_ints (our analysis pipeline)
# rule polyidus has a wildcard contstraint so it can't be used for our pipeline
# so make it higher in the priority

ruleorder: vifi > polyidus > verse 

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
	outpath = row['outdir']
	exp = row['experiment']
	dset = row['analysis_condition']
	host = row['host']
	virus = row['virus']
		
	samples = sample_dict[exp]
		
	for samp in samples:
		other_tool_targets.add(
			f"{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.txt"
		)

rule all:
	input:
		set(analysis_summaries),
		other_tool_targets
	
include: "snakemake_rules/trim.smk"

include: "snakemake_rules/polyidus.smk"

include: "snakemake_rules/seeksv.smk"

include: "snakemake_rules/verse.smk"

include: "snakemake_rules/vifi.smk"
