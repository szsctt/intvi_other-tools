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
		cp $(realpath {input.genome}) $(realpath {output.genome})
		"""
		
#rule vseq_toolkit_config:
#	input:
#		combined_idx = rules.host_virus_index_seeksv.output.idx,
#		vec_idx = multiext("{outpath}/refs/bwa/{virus}/{virus}.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
#		fastq1 = lambda wildcards: get_input_reads(wildcards, analysis_df, 1),
#		fastq2 = lambda wildcards: get_input_reads(wildcards, analysis_df, 2),
#	output:
#		config = "{outpath}/{dset}/vseq_toolkit/{samp}.{host}.{virus}.config.txt"
#	params:
#		adapter1 = lambda wildcards: analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'adapter_1'),
#		adapter2 = lambda wildcards: analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'adapter_2'),	
#		vec_idx = lambda wildcards, input: os.path.abspath(os.path.splitext(input.vec_idx[0])[0]),
#		combined_idx = lambda wildcards, input: os.path.splitext(input.combined_idx[0])[0],
#		qual = lambda wildcards: int(analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'qual')),
#		lenPer = lambda wildcards: int(analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'lenPer')),
#		mode = lambda wildcards: analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'mode'),
#		vecVecFusion = lambda wildcards: analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'vecVecFusion'),
#		stringencyVec = lambda wildcards: analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'stringencyVec'),
#		UMthresholdVec = lambda wildcards: analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'UMthresholdVec'),
#		minMapSpanVec = lambda wildcards: int(analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'minMapSpanVec')),
#		distVecVec = lambda wildcards: int(analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'distVecVec')),
#		opVecVec = lambda wildcards: int(analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'opVecVec')),
#		idenVecVec = lambda wildcards: int(analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'idenVecVec')),
#		stringencyVecGen = lambda wildcards: analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'stringencyVecGen'),
#		UMthresholdVecGen = lambda wildcards: analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'UMthresholdVecGen'),
#		minMapSpanVecGen = lambda wildcards: int(analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'minMapSpanVecGen')),
#		distVecGen = lambda wildcards: int(analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'distVecGen')),
#		opVecGen = lambda wildcards: int(analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'opVecGen')),
#		idenVecGen = lambda wildcards: int(analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'idenVecGen')),
#		clusterRange = lambda wildcards: int(analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'clusterRange')),
#		annoTable = lambda wildcards: analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'host_table')
#	container:
#		"docker://szsctt/vseq:1"
#	shell:
#		"""
#		# General Parameters: Always set these parameters for executing any of the three MODES.
#		echo "file1= $(realpath {input.fastq1})" > {output.config}
#		echo "file2= $(realpath {input.fastq2})" >> {output.config}
#		echo "outDir= $(realpath {wildcards.outpath}/{wildcards.dset}/vseq_toolkit/{wildcards.samp}.{wildcards.host}.{wildcards.virus}/)/" >> {output.config}
#		echo "bin= $VSeqToolkit/scripts/" >> {output.config}
#		echo "qua={params.qual}" >> {output.config}
#		echo "lenPer={params.lenPer}" >> {output.config}
#		echo "adapter1={params.adapter1}" >> {output.config}
#		echo "adapter2={params.adapter2}" >> {output.config}
#		echo "trimmer= $VSeqToolkit/thirdPartyTools/skewer" >> {output.config}
#		echo "aligner= $VSeqToolkit/thirdPartyTools/bwa" >> {output.config}
#		echo "samtools= $VSeqToolkit/thirdPartyTools/samtools" >> {output.config}
#		echo "mode={params.mode}" >> {output.config}
#		
#		# MODE 1: Contaminants Distribution Analysis
#		echo "contAna=false" >> {output.config}
#		echo "combinedRef= $VSeqToolkit/testDir/testReferenceIndex/referenceTestCombined.fa" >> {output.config}
#		echo "stringencyCont=low" >> {output.config}
#		echo "UMthresholdCont=0.95" >> {output.config}
#		
#		#  MODE 2: Vector-Vector Fusion Analysis
#		echo "vecVecFusion={params.vecVecFusion}" >> {output.config}
#		echo "vecRef= {params.vec_idx}" >> {output.config}
#		echo "stringencyVec={params.stringencyVec}" >> {output.config}
#		echo "UMthresholdVec={params.UMthresholdVec}" >> {output.config}
#		echo "minMapSpanVec={params.minMapSpanVec}" >> {output.config}
#		echo "distVecVec={params.distVecVec}" >> {output.config}
#		echo "opVecVec={params.opVecVec}" >> {output.config}
#		echo "idenVecVec={params.idenVecVec}" >> {output.config}
#		
#		# MODE 3: Vector-Host Fusion Analysis 
#		echo "vecGenIS=true" >> {output.config}
#		echo "vecRef= {params.vec_idx}" >> {output.config}
#		echo "vecGenRef= $(realpath {params.combined_idx})" >> {output.config}
#		echo "stringencyVecGen=high" >> {output.config}
#		echo "UMthresholdVecGen={params.UMthresholdVecGen}" >> {output.config}
#		echo "minMapSpanVecGen={params.minMapSpanVecGen}" >> {output.config}
#		echo "distVecGen={params.distVecGen}" >> {output.config}
#		echo "opVecGen={params.opVecGen}" >> {output.config}
#		echo "idenVecGen={params.idenVecGen}" >> {output.config}
#		echo "clusterRange={params.clusterRange}" >> {output.config}
#		echo "annoTable= $VSeqToolkit/testDir/testReferenceIndex/refSeqUCSCTablehg38.txt" >> {output.config}
#		echo "bedtools= $VSeqToolkit/thirdPartyTools/bedtools" >> {output.config}
#		"""
		
		
rule vseq_toolkit_config_template:
	input:
		combined_idx = rules.host_virus_index_seeksv.output.idx,
		vec_idx = multiext("{outpath}/refs/bwa/{virus}/{virus}.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
		fastq1 = lambda wildcards: get_input_reads(wildcards, analysis_df, 1),
		fastq2 = lambda wildcards: get_input_reads(wildcards, analysis_df, 2),
	output:
		config = "{outpath}/{dset}/vseq_toolkit/{samp}.{host}.{virus}.config.txt"
	params:
		fastq1 = lambda wilcards, input: os.path.abspath(input.fastq1),
		fastq2 = lambda wilcards, input: os.path.abspath(input.fastq2),
		out_dir = lambda wildcards: os.path.abspath(f"{wildcards.outpath}/{wildcards.dset}/vseq_toolkit/{wildcards.samp}.{wildcards.host}.{wildcards.virus}"),
		adapter1 = lambda wildcards: analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'adapter_1'),
		adapter2 = lambda wildcards: analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'adapter_2'),	
		vec_idx = lambda wildcards, input: os.path.abspath(os.path.splitext(input.vec_idx[0])[0]),
		combined_idx = lambda wildcards, input: os.path.abspath(os.path.splitext(input.combined_idx[0])[0]),
		qual = lambda wildcards: int(analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'qual')),
		lenPer = lambda wildcards: int(analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'lenPer')),
		mode = lambda wildcards: analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'mode'),
		vecVecFusion = lambda wildcards: analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'vecVecFusion'),
		stringencyVec = lambda wildcards: analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'stringencyVec'),
		UMthresholdVec = lambda wildcards: analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'UMthresholdVec'),
		minMapSpanVec = lambda wildcards: int(analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'minMapSpanVec')),
		distVecVec = lambda wildcards: int(analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'distVecVec')),
		opVecVec = lambda wildcards: int(analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'opVecVec')),
		idenVecVec = lambda wildcards: int(analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'idenVecVec')),
		stringencyVecGen = lambda wildcards: analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'stringencyVecGen'),
		UMthresholdVecGen = lambda wildcards: analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'UMthresholdVecGen'),
		minMapSpanVecGen = lambda wildcards: int(analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'minMapSpanVecGen')),
		distVecGen = lambda wildcards: int(analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'distVecGen')),
		opVecGen = lambda wildcards: int(analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'opVecGen')),
		idenVecGen = lambda wildcards: int(analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'idenVecGen')),
		clusterRange = lambda wildcards: int(analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'clusterRange')),
		annoTable = lambda wildcards: analysis_df_tool_value(wildcards, analysis_df, 'vseq_toolkit', 'host_table')
	container:
		"docker://szsctt/vseq:1"
	shell:
		"""
		cp /var/work/VSeq-Toolkit/config.test.txt {output.config}
		sed -i 's#file1= $VSeqToolkit/testDir/testData/testDataCombined.R1.fastq.gz#file1= {params.fastq1}#g' {output.config}
		sed -i 's#file2= $VSeqToolkit/testDir/testData/testDataCombined.R2.fastq.gz#file2= {params.fastq2}#g' {output.config}
		sed -i 's#outDir= $VSeqToolkit/testDir/testResultsCheck/#outDir= {params.out_dir}/#g' {output.config}
		sed -i 's/qua=20/qua={params.qual}/g' {output.config}
		sed -i 's/lenPer=50/lenPer={params.lenPer}/g' {output.config}
		sed -i 's/adapter1=GATCGGAAGAGCACACGTCTGAACTCCAGTCAC/adapter1={params.adapter1}/g' {output.config}
		sed -i 's/adapter2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT/adapter2={params.adapter2}/g' {output.config}
		sed -i 's/mode=default/mode={params.mode}/g' {output.config}
		sed -i 's/contAna=true/contAna=false/g' {output.config}		
		sed -i 's/vecVecFusion=true/vecVecFusion={params.vecVecFusion}/g' {output.config}
		sed -i 's#vecRef= $VSeqToolkit/testDir/testReferenceIndex/vector1.fa#vecRef= {params.vec_idx}#g' {output.config}
		sed -i 's/stringencyVec=low/stringencyVec={params.stringencyVec}/g'  {output.config}
		sed -i 's/UMthresholdVec=0.95/UMthresholdVec={params.UMthresholdVec}/g' {output.config}
		sed -i 's/minMapSpanVec=20/minMapSpanVec={params.minMapSpanVec}/g' {output.config}		
		sed -i 's/distVecVec=10/distVecVec={params.distVecVec}/g' {output.config}		
		sed -i 's/opVecVec=5/opVecVec={params.opVecVec}/g' {output.config}		
		sed -i 's/idenVecVec=95/idenVecVec={params.idenVecGen}/g' {output.config}		
		sed -i 's#vecGenRef= $VSeqToolkit/testDir/testReferenceIndex/hg38chr22Vector1.fa#vecGenRef= {params.combined_idx}#g' {output.config}
		sed -i 's/stringencyVecGen=low/stringencyVecGen={params.stringencyVecGen}/g' {output.config}
		sed -i 's/UMthresholdVecGen=0.95/UMthresholdVecGen={params.UMthresholdVecGen}/g' {output.config}
		sed -i 's/minMapSpanVecGen=20/minMapSpanVecGen={params.minMapSpanVecGen}/g' {output.config}
		sed -i 's/distVecGen=10/distVecGen={params.distVecGen}/g' {output.config}
		sed -i 's/opVecGen=/opVecGen={params.opVecGen}/g' {output.config}
		sed -i 's/idenVecGen=/idenVecGen={params.idenVecGen}/g' {output.config}
		sed -i 's/clusterRange=3/clusterRange={params.clusterRange}/g' {output.config}
		sed -i 's#annoTable= $VSeqToolkit/testDir/testReferenceIndex/refSeqUCSCTablehg38.txt#annoTable= {params.annoTable}#g' {output.config}
		"""
		
rule vseq_toolkit:
	input:
		combined_idx = rules.host_virus_index_seeksv.output.idx,
		combined_genome = rules.host_virus_index_seeksv.output.fa,
		vec_idx = multiext("{outpath}/refs/bwa/{virus}/{virus}.fa", ".amb", ".ann", ".bwt", ".pac", ".sa"),
		vec_genome = "{outpath}/refs/bwa/{virus}/{virus}.fa",
		fastq1 = lambda wildcards: get_input_reads(wildcards, analysis_df, 1),
		fastq2 = lambda wildcards: get_input_reads(wildcards, analysis_df, 2),
		config = rules.vseq_toolkit_config_template.output.config
	output:
		csv = "{outpath}/{dset}/vseq_toolkit/{samp}.{host}.{virus}/ISGenomeVector.csv"
	container:
		"docker://szsctt/vseq:1"
	threads: 10 # hard-coded into scripts that run bwa-mem?
	resources:
		mem_mb= lambda wildcards, attempt, input: resources_list_with_min_and_max((input.combined_genome, input.vec_genome), attempt, 20000),
		time = lambda wildcards, attempt: ('2:00:00', '24:00:00', '24:00:00', '7-00:00:00')[attempt - 1],
		nodes = 1	
	shell:
		"""
		#CONFIG=$(realpath {input.config})
		#cd $(dirname {output.csv})
		
		#perl /var/work/VSeq-Toolkit/scripts/VSeq-TK.pl -c $CONFIG
		
		#if [ ! -e $(basename {output.csv}) ]; then
		#	touch $(basename {output.csv})
		#fi
		
		perl /var/work/VSeq-Toolkit/scripts/VSeq-TK.pl -c {input.config}
		
		if [ ! -e {output.csv} ]; then
			touch {output.csv}
		fi
		
		"""		
		
				
