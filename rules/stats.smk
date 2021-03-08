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


