
rule mvBamBai:
  input:
    bam = 'bbmap/{sample}.bam'
    #bam = lambda wildcards: os.path.join('bbmap', "{ids}").format(ids=wildcards.id)
    #bam = lambda wildcards: glob('bbmap/*.bam') 
    #bam = lambda wildcards: ids[ids['sample'] == wildcards.idNest]['sample'].tolist(),
    #bai = lambda wildcards: ids[ids['sample'] == wildcards.id]['sample'].tolist(),
    ##bam = lambda wildcards: getBamHome(wildcards.idNest, 'bam'),
    #bai = lambda wildcards: getBamHome(wildcards.id, 'bam.bai')
    #bam = lambda wildcards: glob('bbmap/{sample}.bam'),
    #bai = 'bbmap/{wildcards.sample}.bam.bai'
  output:
    #bam = 'bbmap/all/{iddict[sample][0]}.bam',
    #bai = 'bbmap/all/{idNest}.bam.bai'
  #run:
  #  renameBamBai('bbmap', iddict)
  message: """--- Renaming bam and bai files to <id_nest>.bam,bam.bai from sample info ---"""
  shell:
    """
    echo {input.bam} &
    echo {output.bam}
    """

#echo {input.bam}.bam \n {input.bam}.bai


#rule callVars1:
#  input:
#    #sam = lambda wildcards: getAllSams(wildcards.allTAs),
#    ref = config['ref'],
#    id_list = ""
#  output:
#    vcf = "vars/ants_initial.vcf"
#  threads: 24
#  message: """--- Calling variants (bbtools) for T. alpestre samples ---"""
#  shell: 
#    """
#    callvariants.sh t={threads} list={input.id_list ref={input.ref} ploidy=2 multisample out={output.vcf}
#    """
#
#
#rule qualCalc:
#  input:
#    bam = "bbmap/{sample}_trm.bam",
#    vcf = "vars/ants_initial.vcf"
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

