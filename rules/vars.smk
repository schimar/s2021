
rule mvBamBai:
  input:
    bam = 'bbmap/{sample}.bam',
    bai = 'bbmap/{sample}.bam.bai'
  output:
    touch('bbmap/{sample}.done')
  message: """--- Renaming bam and bai files to <id_nest>.bam,bam.bai from sample info ---"""
  run:
    renameBamBai({input.bam}, {input.bai})


rule lsBam:
  input:
    #bam = 'bbmap/all/{sample}.bam'
  output:
    #touch('bbmap/all/{sample}.ls.done'),
    'vars/bam.list'
  message: """--- creating sample list for variant calling ---"""
  shell:
    """
    ls bbmap/all/*.bam > vars/bam.list  
    """


rule callVars:
  input:
    ref = config['ref'],
    id_list = "vars/bam.list"
  output:
    vcf = "vars/ta_init.vcf"
  threads: 24
  message: """--- Calling variants (bbtools) for T. alpestre samples ---"""
  shell: 
    """
    /apps/uibk/bin/sysconfcpus -n 24 callvariants.sh -Xmx240g -Xms240g t={threads} list={input.id_list} ref={input.ref} ploidy=2 multisample out={output.vcf}
    """

rule filterVars:
  input:
    'vars/ta_init.vcf'
  output:
    'vars/ta.vcf'
  message: """--- Filtering variants for T. alpestre samples ---"""
  shell: 
    """
    script/vcfBBfilt.py {input} > {output}
    """



#rule qualCalc:
#  input:
#    bam = "bbmap/{sample}.bam",
#    vcf = "vars/ta_init.vcf"
#    ref = config['ref']
#  output:
#    "placeholder_for_output"
#  threads: 12
#  message: """--- Calculating true quality with bbtools ---"""
#  shell:
#    """
#    calctruequality.sh t={threads} vcf={input.vcf} ref={input.ref} id= [[a,b,c]] 
#    """
#
#
#rule bbdukQualRecal:
#  input:
#    bam = 'bbmap/{sample}.bam'
#  output:
#    'what will it be? perhaps simply the same bam files?'
#  threads: 12
#  message: """--- Recalibrating quality of bam files ---"""
#  shell:
#    """
#    
#    """

rule lsBamRecal:
  input: 
  output:
    'vars/bam.recal.list'
  message: """--- creating recalibrated sample list for variant calling ---"""
  shell:
    """
    ls bbmap/all/*.recal.bam > vars/bam.recal.list  
    """

#
#rule callVars2:
#  input:
#  ref = config['ref'],
#    id_list = ""
#  output:
#    vcf = "vars/ants_initial.vcf"
#  threads: 24
#  message: """--- Calling variants (bbtools) for T. alpestre samples ---"""
#  shell: 
#    """
#    callvariants.sh t={threads} list={input.id_list ref={input.ref} ploidy=2 multisample out={output.vcf}
#    """

