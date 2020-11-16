
#####################################################
##################### polyidus ######################
#####################################################

rule bwt2_index:
	input:
		fasta = lambda wildcards: ref_names[wildcards.genome]
	output:
		multiext("{outpath}/references/{genome}/{genome}", 
							".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"
						)
	container:
		"docker://szsctt/polyidus:2"
	params:
		prefix = lambda wildcards, output: os.path.splitext(os.path.splitext(output[0])[0])[0]
	resources:
		mem_mb= lambda wildcards, attempt, input: int(attempt * 5 * (os.stat(input.fasta).st_size/1e6)),
		time = lambda wildcards, attempt: ('2:00:00', '24:00:00', '24:00:00', '7-00:00:00')[attempt - 1],
		nodes = 1
	shell:
		"bowtie2-build {input} {params.prefix}"

rule trim:
	input:
		r1 = lambda wildcards: get_input_reads(wildcards, 1),
		r2 = lambda wildcards: get_input_reads(wildcards, 2)
	output:
		proc_r1 = temp("{outpath}/{dset}/trimmed_reads/{samp}.1.fastq.gz"),
		proc_r2 = temp("{outpath}/{dset}/trimmed_reads/{samp}.2.fastq.gz")
	conda:	
		"../envs/seqprep.yml"
	container:
		"docker://szsctt/seqprep:1"
	params:
		A = lambda wildcards: analysis_df_value(wildcards, 'adapter_1'),
		B = lambda wildcards: analysis_df_value(wildcards, 'adapter_2')
	shell:
		"""
		SeqPrep -A {params.A} -B {params.B} -f {input.r1} -r {input.r2} -1 {output.proc_r1} -2 {output.proc_r2}
		"""


rule polyidus:
	input:
		fastq1 = lambda wildcards: get_polyidus_reads(wildcards, 1),
		fastq2 = lambda wildcards: get_polyidus_reads(wildcards, 2),
		host_idx = multiext("{outpath}/references/{host}/{host}", 
							".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"
						),
		virus_idx =  multiext("{outpath}/references/{virus}/{virus}", 
							".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"
						)
	output:
		temp_ints = temp("{outpath}/{dset}/polyidus.{host}.{virus}.{samp}/results/exactHpvIntegrations.tsv"),
		ints = "{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.txt",
		fake_merged = temp("{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.merged.bed")
	params:
		output = lambda wildcards, output: os.path.dirname(os.path.dirname(output.temp_ints)),
		host_idx = lambda wildcards, input: os.path.splitext(os.path.splitext(input.host_idx[0])[0])[0],
		virus_idx = lambda wildcards, input: os.path.splitext(os.path.splitext(input.virus_idx[0])[0])[0],
	resources:
		mem_mb= lambda wildcards, attempt: int(attempt * 10000),
		time = lambda wildcards, attempt: ('2:00:00', '24:00:00', '24:00:00', '7-00:00:00')[attempt - 1]
	container:
		"docker://szsctt/polyidus:2"
	wildcard_constraints:
		dset = ".+_polyidus\d+"
	shell:
		"""
		rm -rf {params.output}/*
		
		python3 /usr/src/app/src/polyidus.py \
		{params.host_idx} \
		{params.virus_idx} \
		--fastq {input.fastq1} {input.fastq2} \
		--outdir {params.output}
		
		cp {output.temp_ints} {output.ints}
		cp {output.temp_ints} {output.fake_merged}
		"""
		
		
def analysis_df_value(wildcards, column_name):
	
	# get a value from the row of the df corresponding to this analysis condition
	unique = f"{wildcards.dset}"

	return analysis_df.loc[(analysis_df['analysis_condition'] == unique).idxmax(), column_name] 


def get_input_reads(wildcards, read_num):
	assert read_num in (1, 2)
	if read_num == 1:
		return f"{analysis_df_value(wildcards, 'read_folder')}/{wildcards.samp}{analysis_df_value(wildcards, 'R1_suffix')}"
	if read_num == 2:
		return f"{analysis_df_value(wildcards, 'read_folder')}/{wildcards.samp}{analysis_df_value(wildcards, 'R2_suffix')}"	


def get_polyidus_reads(wildcards, read_num):

	if analysis_df_value(wildcards, 'trim') == 1:
		if read_num == 1:
			return rules.trim.output.proc_r1
		else:
			return rules.trim.output.proc_r2
	else:
		if read_num == 1:
			return get_input_reads(wildcards, 1)
		else:
			return get_input_reads(wildcards, 2)
