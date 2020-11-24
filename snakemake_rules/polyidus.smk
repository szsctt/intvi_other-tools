
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
		mem_mb= lambda wildcards, attempt, input: resources_list_with_min_and_max((input.fasta, ), attempt),
		time = lambda wildcards, attempt: ('2:00:00', '24:00:00', '24:00:00', '7-00:00:00')[attempt - 1],
		nodes = 1
	shell:
		"bowtie2-build {input} {params.prefix}"

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
		mem_mb= lambda wildcards, attempt, input: resources_list_with_min_and_max(input.host_idx, attempt),
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
		
		

