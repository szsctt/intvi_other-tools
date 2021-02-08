#####################################################
################### VSeq-Toolkit ####################
#####################################################

rule bwa_index:
	input:
		genome = lambda wildcards: ref_names[wildcards.genome],
	output:
		idx = multiext("{outpath}/refs/bwa/{genome}/{genome}.fa", ".amb", ".ann", ".bwt", ".pac", ".sa")
	params:
		prefix = lambda wildcards, output: os.path.splitext(output[0])[0]
	container:
		"docker://szsctt/vseq:1"
	resources:
		mem_mb= lambda wildcards, attempt, input: resources_list_with_min_and_max((input.genome, ), attempt),
		time = lambda wildcards, attempt: ('2:00:00', '24:00:00', '24:00:00', '7-00:00:00')[attempt - 1],
		nodes = 1	
	shell:
		"""
		 $VSeqToolkit/thirdPartyTools/bwa index -p {params.prefix} {input.genome}
		"""
rule link_virus:
	input:
		genome = lambda wildcards: ref_names[wildcards.genome],
	output:
		genome = "{outpath}/refs/bwa/{genome}/{genome}.fa"
	container:
		"docker://szsctt/vseq:1"
	shell:
		"""
		ln -sf $(realpath {input.genome}) $(realpath {output.genome})
		"""
		
rule vseq_toolkit_config:
	input:
		combined_idx = rules.host_virus_index_seeksv.output.idx,
		vec_idx = multiext("{outpath}/refs/bwa/{virus}/{virus}.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
		fastq1 = lambda wildcards: get_input_reads(wildcards, analysis_df, 1),
		fastq2 = lambda wildcards: get_input_reads(wildcards, analysis_df, 2),
	output:
		config = "{outpath}/{dset}/vseq_toolkit/{samp}.{host}.{virus}/config.txt"
	params:
		adapter1 = lambda wildcards: analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'adapter_1'),
		adapter2 = lambda wildcards: analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'adapter_2'),
		vec_idx = lambda wildcards, input: os.path.splitext(input.vec_idx[0])[0],
		combined_idx = lambda wildcards, input: os.path.splitext(input.combined_idx[0])[0],
	container:
		"docker://szsctt/vseq:1"
	shell:
		"""
		# General Parameters: Always set these parameters for executing any of the three MODES.
		echo "file1= $(realpath {input.fastq1})" > {output.config}
		echo "file2= $(realpath {input.fastq2})" >> {output.config}
		echo "outDir= $(realpath {wildcards.outpath}/{wildcards.dset}/vseq_toolkit/{wildcards.samp}.{wildcards.host}.{wildcards.virus}/)/" >> {output.config}
		echo "bin= $VSeqToolkit/scripts/" >> {output.config}
		echo "qua=20" >> {output.config}
		echo "lenPer=50" >> {output.config}
		echo "adapter1={params.adapter2}" >> {output.config}
		echo "adapter2={params.adapter2}" >> {output.config}
		echo "trimmer= $VSeqToolkit/thirdPartyTools/skewer" >> {output.config}
		echo "aligner= $VSeqToolkit/thirdPartyTools/bwa" >> {output.config}
		echo "samtools= $VSeqToolkit/thirdPartyTools/samtools" >> {output.config}
		echo "mode=default" >> {output.config}
		
		# MODE 1: Contaminants Distribution Analysis
		echo "contAna=false" >> {output.config}
		echo "combinedRef= $VSeqToolkit/testDir/testReferenceIndex/referenceTestCombined.fa" >> {output.config}
		echo "stringencyCont=low" >> {output.config}
		echo "UMthresholdCont=0.95" >> {output.config}
		
		#  MODE 2: Vector-Vector Fusion Analysis
		echo "vecVecFusion=false" >> {output.config}
		echo "vecRef= $(realpath {params.vec_idx})" >> {output.config}
		echo "stringencyVec=low" >> {output.config}
		echo "UMthresholdVec=0.95" >> {output.config}
		echo "minMapSpanVec=20" >> {output.config}
		echo "distVecVec=10" >> {output.config}
		echo "opVecVec=5" >> {output.config}
		echo "idenVecVec=95" >> {output.config}
		
		# MODE 3: Vector-Host Fusion Analysis 
		echo "vecGenIS=true" >> {output.config}
		echo "vecRef= $(realpath {params.vec_idx})" >> {output.config}
		echo "vecGenRef= $(realpath {params.combined_idx})" >> {output.config}
		echo "stringencyVecGen=high" >> {output.config}
		echo "UMthresholdVecGen=0.95" >> {output.config}
		echo "minMapSpanVecGen=20" >> {output.config}
		echo "distVecGen=10" >> {output.config}
		echo "opVecGen=5" >> {output.config}
		echo "idenVecGen=95" >> {output.config}
		echo "clusterRange=3" >> {output.config}
		echo "annoTable= $VSeqToolkit/testDir/testReferenceIndex/refSeqUCSCTablehg38.txt" >> {output.config}
		echo "bedtools= $VSeqToolkit/thirdPartyTools/bedtools" >> {output.config}
		"""
		
rule vseq_toolkit:
	input:
		combined_idx = rules.host_virus_index_seeksv.output.idx,
		combined_genome = rules.host_virus_index_seeksv.output.fa,
		vec_idx = multiext("{outpath}/refs/bwa/{virus}/{virus}.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
		vec_genome = "{outpath}/refs/bwa/{virus}/{virus}.fa",
		fastq1 = lambda wildcards: get_input_reads(wildcards, analysis_df, 1),
		fastq2 = lambda wildcards: get_input_reads(wildcards, analysis_df, 2),
		config = rules.vseq_toolkit_config.output.config
	output:
		csv = "{outpath}/{dset}/vseq_toolkit/{samp}.{host}.{virus}/ISGenomeVector.csv"
	container:
		"docker://szsctt/vseq:1"
	shell:
		"""
		perl $VSeqToolkit/scripts/VSeq-TK.pl -c {input.config}
		"""		
		
				
