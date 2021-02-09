include: "rules/common.smk"



# -----------------------------------------------


rule all:
    input:
        expand('raw/qc/fastqc/{sample}_{pair}_001_fastqc.{ext}', sample=sample_names, pair=['R1', 'R2'], ext=['html', 'zip']),
        expand('trm/{sample}_{pair}_trm.fq.gz', sample=sample_names, pair=['R1', 'R2']), 
	    expand('bbmap/{sample}_trm.bam', sample=sample_names),
	    expand('bbmap/{sample}_trm.bam.bai', sample=sample_names),
		'ref/genome/1/summary.txt'



# -----------------------------------------------


include: "rules/hts.smk"
include: "rules/vars.smk"




