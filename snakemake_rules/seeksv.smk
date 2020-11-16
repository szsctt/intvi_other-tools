
rule host_virus_index_seeksv:
	input:
		virus = lambda wildcards: ref_names[wildcards.virus],
		host = lambda wildcards: ref_names[wildcards.host]
	output:
		fa = "{outpath}/seeksv_refs/data/{virus}/{host}_{virus}.fas",
		idx = multiext("{outpath}/seeksv_refs/data/{virus}/{host}_{virus}.fas", ".amb", ".ann", ".bwt", ".pac", ".sa")
	container:
		"docker://szsctt/seeksv:1"
	resources:
		mem_mb= lambda wildcards, attempt, input: int(attempt * 5 * (os.stat(input.host).st_size/1e6 + os.stat(input.virus).st_size/1e6)),
		time = lambda wildcards, attempt: ('2:00:00', '24:00:00', '24:00:00', '7-00:00:00')[attempt - 1],
		nodes = 1
	shell:
		"""
		cat {input.virus} {input.host} > {output.fa}
		bwa index {output.fa}
		"""
		
rule align_seeksv_all:
	input:
		fastq1 = lambda wildcards: get_polyidus_reads(wildcards, 1),
		fastq2 = lambda wildcards: get_polyidus_reads(wildcards, 2),
		idx = rules.host_virus_index_seeksv.output.idx
	output:
		bam = "{outpath}/{dset}/aln/{samp}.{host}.{virus}.bam",
		bai = "{outpath}/{dset}/aln/{samp}.{host}.{virus}.bam.bai"
	params:
		prefix = lambda wildcards, input: os.path.splitext(input.idx[0])[0]	
	threads: 8
	resources:
		mem_mb = 10000
	container:
		"docker://szsctt/seeksv:1"	
	shell:
		"""
		bwa mem {params.prefix} {input.fastq1} {input.fastq2} | samtools sort -o {output.bam} -
		samtools index {output.bam}
		"""
		
rule dedup_seeksv:
	input:
		rules.align_seeksv_all.output.bam
	output:
		bam = "{outpath}/{dset}/aln/{samp}.{host}.{virus}.dup.bam",
		metrics = "{outpath}/{dset}/aln/{samp}.{host}.{virus}.dup.txt"
	params:
		mem_gb_sort = lambda wildcards, resources: int(resources.mem_mb / 1e3 / 3),
		mem_gb_dup = lambda wildcards, resources: int(resources.mem_mb / 1e3 / 3),
	resources:
		mem_mb = 10000
	container:
		"docker://szsctt/seeksv:1"	
	shell:
		"""
		java -Xmx{params.mem_gb_sort}g -jar ${{PICARD}} SortSam \
			I=/dev/stdin \
			VALIDATION_STRINGENCY=LENIENT \
			COMPRESSION_LEVEL=0 \
			O=/dev/stdout \
			SORT_ORDER=queryname | java -Xmx{params.mem_gb_dup}g -jar ${{PICARD}} MarkDuplicates \
			I=/dev/stdin \
			VALIDATION_STRINGENCY=LENIENT \
			METRICS_FILE={output.metrics} \
			O={output.bam}
		"""
		
def get_seeksv_alignment(wildcards):
	# if we want to do dedupliation
	if analysis_df_value(wildcards, 'dedup') == 1:
		return rules.dedup_seeksv.output.bam
	
	# no deduplication
	return rules.align_seeksv_all.output.bam
		
rule sort_seeksv:
	input:
		bam = lambda wildcards: get_seeksv_alignment(wildcards)
	output:
		bam = "{outpath}/{dset}/aln/{samp}.{host}.{virus}.sorted.bam",
		bai = "{outpath}/{dset}/aln/{samp}.{host}.{virus}.sorted.bam.bai"
	resources:
		mem_mb = 10000
	params:
		mem_gb = lambda wildcards, resources: int(resources.mem_mb / 1e3) - 1
	container:
		"docker://szsctt/seeksv:1"	
	shell:
		"""
		java -Xmx{params.mem_gb}g -jar ${{PICARD}} SortSam \
			I={input.bam} \
			VALIDATION_STRINGENCY=LENIENT \
			COMPRESSION_LEVEL=0 \
			O={output.bam} \
			SORT_ORDER=coordinate
			
		samtools index {output.bam}
		"""
		
rule seeksv_getclip:
	input:
		bam = rules.sort_seeksv.output.bam,
		bai = rules.sort_seeksv.output.bai
	output: 
		clip_fq = "{outpath}/{dset}/clipped_reads/{samp}.{host}.{virus}.clip.fq.gz",
		clip = "{outpath}/{dset}/clipped_reads/{samp}.{host}.{virus}.clip.gz",
		unmapped_1 = "{outpath}/{dset}/clipped_reads/{samp}.{host}.{virus}.unmapped_1.fq.gz",
		unmapped_2 = "{outpath}/{dset}/clipped_reads/{samp}.{host}.{virus}.unmapped_2.fq.gz",	
	params:
		prefix = lambda wildcards, output: os.path.splitext(os.path.splitext(output.clip)[0])[0]
	container:
		"docker://szsctt/seeksv:1"
	shell:
		"""
		/var/work/seeksv/seeksv getclip -o {params.prefix} {input.bam}
		"""
		
rule align_seeksv_clip:
	input:
		fq = rules.seeksv_getclip.output.clip_fq,
		idx = rules.host_virus_index_seeksv.output.idx
	output:
		bam = "{outpath}/{dset}/aln/{samp}.{host}.{virus}.clip.bam"
	params:
		prefix = lambda wildcards, input: os.path.splitext(input.idx[0])[0]
	container:
		"docker://szsctt/seeksv:1"
	shell:
		"""
		bwa mem  {params.prefix} {input.fq} | samtools view  -Sb -o {output.bam} -
		"""

rule seeksv:
	input:
		clip = rules.seeksv_getclip.output.clip,
		bam = rules.align_seeksv_all.output.bam,
		bam_clip = rules.align_seeksv_clip.output.bam
	output:
		ints = "{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.txt",
		unmapped = "{outpath}/{dset}/ints/{samp}.{host}.{virus}.unmapped.clip.fq.gz",
		fake_merged = temp("{outpath}/{dset}/ints/{samp}.{host}.{virus}.integrations.merged.bed")
	container:
		"docker://szsctt/seeksv:1"
	wildcard_constraints:
		dset = ".+_seeksv\d+"
	shell:
		"""
		/var/work/seeksv/seeksv getsv {input.bam_clip} {input.bam} {input.clip} {output.ints} {output.unmapped}
		cp {output.ints} {output.fake_merged}
		"""		
