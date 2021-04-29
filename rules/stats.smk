rule vcf2zarr:
  input:
    vcf = 'vars/ta{vartype}Bypop.vcf'
  output:
    'vars/ta{vartype}.zarr/.zgroup'
  message: """--- Converting vcf into zarr format ---"""
  shell:
    """
      python script/vcf2zarr.py {input.vcf}
    """

##echo '{input.vcf}.zarr'


rule alStats:
  input:
    zarr = 'vars/ta{vartype}.zarr/'
  output:
    touch("vars/ta{vartype}/al.done")
  message:
    """--- Creating scikit-allel statistics ---"""
  shell:
    """
    script/al.py {input.zarr}
    """

rule gemma_mg:
  input:
    vcf = 'vars/taSubInDel.{sets}.vcf'
    #vcf = 'vars/ta{vartype}Bypop.vcf'
  output:
    mg = 'vars/taSubInDel/stats/gemma/taSubInDel.{sets}.mg'
    #mg = 'vars/ta{vartype}/stats/gemma/{vartype}.mg'
  message:
    """--- Converting vcf to bimbam format ---"""
  shell:
    """
    script/vcf2mg.py {input.vcf} > {output.mg}
    """

rule gemma_relmat_stdzd:
  input:
    mg = 'vars/taSubInDel/stats/gemma/taSubInDel.{sets}.mg',
    pheno = 'vars/taSubInDel/stats/gemma/taSubInDel.pheno'
    #mg = 'vars/ta{vartype}/stats/gemma/ta{vartype}.mg',
    #pheno = 'vars/ta{vartype}/stats/gemma/ta{vartype}.pheno'
  output:
    relM = 'vars/taSubInDel/stats/gemma/taSubInDel.{sets}.stdzd.relmat'
    #relM = 'vars/ta{vartype}/stats/gemma/ta{vartype}.stdzd.relmat'
  message:
    """--- Calculating standardized relatedness matrix with gemma ---"""
  shell:
    """
    gemma -g {input.mg} -p {input.pheno} -gk 2 -o {output.relM}
    """

rule gemma_relmat_ctrd:
  # NOTE: see **manual 4.4.2** on details regarding centered vs standardized rel.matrix
  input:
    mg = 'vars/taSubInDel/stats/gemma/taSubInDel.{sets}.mg',
    pheno = 'vars/taSubInDel/stats/gemma/taSubInDel.pheno'
    #mg = 'vars/ta{vartype}/stats/gemma/ta{vartype}.mg',
    #pheno = 'vars/ta{vartype}/stats/gemma/ta{vartype}.pheno'
  output:
    touch('vars/taSubInDel/stats/gemma/relmat/ctrd.done'),
    relM = 'taSubInDel'
        #relM = 'vars/ta{vartype}/stats/gemma/ta{vartype}.ctrd.relmat'
  message:
    """--- Calculating centered relatedness matrix with gemma ---"""
  shell:
    """
    cd vars/taSubInDel/stats/gemma/
    gemma -g {input.mg} -p {input.pheno} -gk 1 -o {output.relM}
    cd ../../../../
    """

rule get_vartype:
  input:
    vcf = 'vars/taSubInDel.{sets}.vcf'
  output:
    txt = 'vars/taSubInDel/stats/gemma/taSubInDel.{sets}.typ.txt'
  message: 
    """--- Getting variant type for vcf ---"""
  shell:
    """
    grep -v '#' {input.vcf} | egrep -o 'TYP=[A-Z]+' | cut -f2 -d$'=' > {output.txt}
    """


rule sub_seg_scafbp:
  input:
    vcf = 'vars/taSubInDelBypop.vcf',
    scafbp = 'vars/taSubInDel/stats/gemma/vars_seg.gemma.scafbp'
  output:
    vcf = 'vars/taSubInDel.seg.vcf'
  message:
    """--- Taking subset of LD-pruned variants ---"""
  shell:
    """
    script/sub_vcf_scafbp.py {input.vcf} {input.scafbp} > {output.vcf}
    """

rule sub_ldp_scafbp:
  input:
    vcf = 'vars/taSubInDelBypop.vcf',
    scafbp = 'vars/taSubInDel/stats/al/pca/ld_prunedVars.scafbp.txt'
  output:
    vcf = 'vars/taSubInDel.ldp.vcf'
  message:
    """--- Taking subset of LD-pruned variants ---"""
  shell:
    """
    script/sub_vcf_scafbp.py {input.vcf} {input.scafbp} > {output.vcf}
    """


rule ngsRelate:
  input:
    vcf = 'vars/taSubInDelBypop.vcf'
  output:
    stats = 'vars/taSubInDel/stats/ngsRelate/stats.txt'
  threads: 12
  message:
    """--- Calculating relatedness (etc.) with ngsRelate2 ---"""
  shell:
    """
    module load ngsrelate/2.0
    ngsrelate -h {input.vcf} -p 12 -T GT -c 1 -O {output.stats} -z vars/ids.txt
    """

rule relStats:
  input:
    zarr = 'vars/taSubInDel.zarr/'
  output:
    touch("vars/taSubInDel/nrel.done")
  message:
    """--- Creating ngsRelate summary statistics & figures ---"""
  shell:
    """
    script/nrel.py {input.zarr}
    """

