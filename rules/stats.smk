rule vcf2zarr:
  input:
    vcf = 'vars/ta{vartype}.vcf'
  output:
    'vars/ta{vartype}.zarr/.zgroup'
  message: """--- Converting vcf into zarr format ---"""
  shell:
    """
    echo '{input.vcf}.zarr'
    script/vcf2zarr.py {input.vcf}
    """


rule allelStats:
  input:
    vcf = 'zarr path??'
  output:
    '{vartype}folder/someFigs{.pdf,eps}'
  message:
    """--- Creating scikit-allel statistics ---"""
  shell:
    """
    allel.py
    """


