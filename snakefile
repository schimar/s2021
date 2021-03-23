include: "rules/common.smk"



# -----------------------------------------------


rule all:
    input:
        #expand('raw/qc/fastqc/{sample}_{pair}_001_fastqc.{ext}', sample=sample_names, pair=['R1', 'R2'], ext=['html', 'zip']),
        #expand('trm/{sample}_{pair}.fq.gz', sample=sample_names, pair=['R1','R2']), 
	    #expand('bbmap/{sample}.{ext}', sample=sample_names, ext=['bam', 'bam.bai']),
	    #'ref/genome/1/summary.txt',
	    #expand('bbmap/{sample}.bam', sample=ids['sample']),
	    #expand('bbmap/{sample}.done', sample=ids['sample']),
	    ##expand('bbmap/{id}.{ext}', id=ids['sample'], ext=['bam', 'bam.bai']),
	    #'vars/bam.list',
	    #'vars/ta_init.vcf',
	    #expand('vars/ta{type}.vcf', type=['SubInDel', 'InDel', 'Sub']),
        expand('vars/ta{type}.zarr/.zgroup', type=['SubInDel', 'InDel', 'Sub']),
        'vars/tapopsOrd.txt',
        expand('vars/ta{vartype}Bypop.vcf', vartype=['SubInDel', 'InDel', 'Sub']),
	    expand('vars/ta{vartype}/al.done', vartype=['SubInDel', 'InDel', 'Sub']),
	    #expand('vars/ta{vartype}/stats/gemma/vars_seg.gemma.scafbp', vartype=['SubInDel', 'InDel', 'Sub']), ### output for rule alStats and input for rule gemma
	    #expand('vars/ta{vartype}/stats/gemma/ta{vartype}.mg', vartype=['SubInDel', 'InDel', 'Sub']),
	    #expand('vars/ta{vartype}/stats/gemma/ta{vartype}.pheno', vartype=['SubInDel', 'InDel', 'Sub']),
	    expand('vars/ta{vartype}/stats/gemma/ta{vartype}.ctrd.relmat', vartype=['SubInDel', 'InDel', 'Sub']),
	    expand('vars/ta{vartype}/stats/gemma/ta{vartype}.stdzd.relmat', vartype=['SubInDel', 'InDel', 'Sub']),







# -----------------------------------------------


include: "rules/hts.smk"
include: "rules/vars.smk"
include: "rules/stats.smk"



