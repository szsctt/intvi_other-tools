
rule host_virus_index_seeksv:
	input:
		virus = lambda wildcards: ref_names[wildcards.virus],
		host = lambda wildcards: ref_names[wildcards.host]
	output:
		fa = "{outpath}/vifi_refs/data/{virus}/{host}_{virus}.fas",
		idx = multiext("{outpath}/vifi_refs/data/{virus}/{host}_{virus}.fas", ".amb", ".ann", ".bwt", ".pac", ".sa")
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
		mem_gb_sort = lambda wildcards, resources: int(resources.mem_mb / 1e3 / 3)
		mem_gb_dup = lambda wildcards, resources: int(resources.mem_mb / 1e3 / 3)
	threads: 8
	resources:
		mem_mb = 10000
	container:
		"docker://szsctt/seeksv:1"	
	shell:
		"""
		bwa mem {input.idx} {input.fastq1} {input.fastq2} |\
		java -Xmx{params.mem_gb_sort}g -jar ${PICARD} SortSam \
			I=/dev/stdin \
			MAX_RECORDS_IN_RAM=7000000 \
			VALIDATION_STRINGENCY=LENIENT \
			COMPRESSION_LEVEL=0 \
			O=/dev/stdout \
			SORT_ORDER=queryname | java -Xmx{mem_gb_dup}g -jar ${PICARD} MarkDuplicates \
			I=/dev/stdin \
			VALIDATION_STRINGENCY=LENIENT \
			O=/dev/stdout | java -Xmx{params.mem_gb_sort}g -jar ${PICARD} SortSam \
			I=/dev/stdin \
			MAX_RECORDS_IN_RAM=7000000 \
			VALIDATION_STRINGENCY=LENIENT \
			COMPRESSION_LEVEL=0 \
			O={output.bam} \
			SORT_ORDER=coordinate
		samtools index {output.bam}
		"""
		
rule seeksv_getclip:
	input:
		bam = "{outpath}/{dset}/aln/{samp}.{host}.{virus}.bam",
		bai = "{outpath}/{dset}/aln/{samp}.{host}.{virus}.bam.bai"
	output: 
		clip_fq = "{outpath}/{dset}/clipped_reads/{samp}.{host}.{virus}.clip.fq.gz",
		clip = "{outpath}/{dset}/clipped_reads/{samp}.{host}.{virus}.clip.gz",
		unmapped_1 = "{outpath}/{dset}/clipped_reads/{samp}.{host}.{virus}.unmapped_1.fq.gz",
		unmapped_2 = "{outpath}/{dset}/clipped_reads/{samp}.{host}.{virus}.unmapped_2.fq.gz",	
	params:
		prefix = lambda wildcards, output: path.splitext(path.splitext(output.clip)[0])[0]
	container:
		"docker://szsctt/seeksv:1"
	shell:
		"""
		seeksv getclip -o {params.prefix} {input.bam}
		"""
		
rule align_seeksv_clip:
	input:
		fq = rules.seeksv_getclip.output.clip_fq,
		idx = rules.host_virus_index_seeksv.output.idx
	output:
		bam = "{outpath}/{dset}/aln/{samp}.{host}.{virus}.clip.bam",
		bai = "{outpath}/{dset}/aln/{samp}.{host}.{virus}.clip.bam.bai"
	params:
		prefix = lambda wildcards, input: path.splitext(input.idx[0])[0]
	container:
		"docker://szsctt/seeksv:1"
	shell:
		"""
		bwa mem  {params.prefix} {input.fq} | samtools view  -Sb -o {output.bam} -
  		samtools index {output.bam}
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
		seeksv getsv {input.bam_clip} {input.bam} {input.clip} {output.ints} {output.unmapped}
		cp {output.ints} {output.fake_merged}
		"""		