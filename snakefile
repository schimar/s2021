include: "rules/common.smk"





######    Rules    #####


rule all:
    input:
        expand('raw/qc/fastqc/{sample}_{pair}_001_fastqc.{ext}', sample=sample_names, pair=['R1', 'R2'], ext=['html', 'zip']),
        expand('trm/{sample}_{pair}_trm.fq.gz', sample=sample_names, pair=['R1', 'R2']), 
	    expand('bbmap/{sample}_trm.bam', sample=sample_names),
	    expand('bbmap/{sample}_trm.bam.bai', sample=sample_names),
		'ref/genome/1/summary.txt'



rule qc:
    input:
        r1 = lambda wildcards: getFqHome(wildcards.sample)[0],
        r2 = lambda wildcards: getFqHome(wildcards.sample)[1]
    output:
        "raw/qc/fastqc/{sample}_R1_001_fastqc.html",
        "raw/qc/fastqc/{sample}_R1_001_fastqc.zip",
        "raw/qc/fastqc/{sample}_R2_001_fastqc.html",
        "raw/qc/fastqc/{sample}_R2_001_fastqc.zip"
    log: "log/{sample}.qc.log.txt"
    #resources:
    #    mem = 1000,
    #    time = 300
    threads: 1
	conda:
	    "envs/s21.yaml"
    message: """--- Quality check of raw data with FastQC before trimming."""
    shell: """
        fastqc -o raw/qc/fastqc/ -f fastq {input.r1} &
        fastqc -o raw/qc/fastqc/ -f fastq {input.r2}
        """


rule trim:
    input:
        r1 = lambda wildcards: getFqHome(wildcards.sample)[0],
        r2 = lambda wildcards: getFqHome(wildcards.sample)[1],
        adapters = config["adapters"]
    output:
        r1trmd = "trm/{sample}_R1_trm.fq.gz",
        r2trmd = "trm/{sample}_R2_trm.fq.gz"
    log: "log/{sample}.bbduk_qual_trim.log.txt"
    #resources:
    #    mem = 1000,
    #    time = 300
    threads: 2
	conda:
	    "envs/s21.yaml"
    message: """--- Quality trimming of fastq files before mapping."""
    shell: 
        """
        bbduk.sh -Xmx1g in1={input.r1} in2={input.r2} out1={output.r1trmd} out2={output.r2trmd} trimq=6 qtrim=r hdist=1 bhist=trm/hist/{wildcards.sample}.bhist qhist=trm/hist/{wildcards.sample}.qhist lhist=trm/hist/{wildcards.sample}.lhist tpe tbo #&> log/{wildcards.sample}_bbduk_trm.log
        """

rule refIndex:
	input:
		ref = config['ref']
	output:
		'ref/genome/1/summary.txt'
	shell:
		"""
		bbmap.sh ref={input.ref}
		"""


rule map:
	input:
		#tr1 = "trm/{sample}_R1_trm.fq.gz",
		#tr2 = "trm/{sample}_R2_trm.fq.gz",
		tr1 = lambda wildcards: getTrmHome(wildcards.sample)[0],
		tr2 = lambda wildcards: getTrmHome(wildcards.sample)[1],
		ref = config["ref"]
	output:
		bam = "bbmap/{sample}_trm.bam"
	#log:
	#	"log/{sample}.bbmap.log.txt"
	threads: 6
	message: """--- Mapping reads to reference genome ---"""
	shell:
		"""
		bbmap.sh -Xmx50g t={threads} ref={input.ref} in1={input.tr1} in2={input.tr2} out={output.bam} minid=0.85 rgid={wildcards.sample} rglb=igaDNA rgsm={wildcards.sample} rgpl=ILLUMINA overwrite=f unpigz=t | samtools view -F 4 -Shu /dev/stdin | samtools sort -l 4 - -o {output.bam} #&> log/{wildcards.sample}.bbmap.log.txt
		"""


rule samIndex:
	input: 
		aln = expand('bbmap/{sample}_trm.bam', sample=sample_names)	#"bbmap/{sample}_trm.bam"
	output:
		idx = "bbmap/{sample}_trm.bam.bai"
	threads: 2
	message: """--- Indexing with samtools ---"""
	shell:
		"""
		samtools index {input.aln}
		"""





