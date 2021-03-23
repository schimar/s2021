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
    zarr = 'vars/ta{vartype}.zarr'
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
    vcf = 'vars/ta{vartype}Bypop.vcf'
  output:
    mg = 'vars/ta{vartype}/stats/gemma/{vartype}.mg'
  message:
    """--- Converting {vartype} vcf to bimbam format ---"""
  shell:
    """
    script/vcf2mg.py {input.vcf} > {output.mg}
    """

rule gemma_relmat_stdzd:
  input:
    mg = 'vars/ta{vartype}/stats/gemma/ta{vartype}.mg'
    pheno = 'vars/ta{vartype}/stats/gemma/ta{vartype}.pheno'
  output:
    relM = 'vars/ta{vartype}/stats/gemma/ta{vartype}.stdzd.relmat'
  message:
    """--- Calculating standardized relatedness matrix with gemma ---"""
  shell:
    """
    gemma -g {input.mg} -p {input.pheno} -gk 2 -o {output.relM}
    """

rule gemma_relmat_ctrd:
  # NOTE: see **manual 4.4.2** on details regarding centered vs standardized rel.matrix
  input:
    mg = 'vars/ta{vartype}/stats/gemma/ta{vartype}.mg'
    pheno = 'vars/ta{vartype}/stats/gemma/ta{vartype}.pheno'
  output:
    relM = 'vars/ta{vartype}/stats/gemma/ta{vartype}.ctrd.relmat'
  message:
    """--- Calculating centered relatedness matrix with gemma ---"""
  shell:
    """
    gemma -g {input.mg} -p {input.pheno} -gk 1 -o {output.relM}
    """




