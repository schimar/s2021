rule vcf2zarr:
  input:
    vcf = 'vars/ta{vartype}.vcf'
  output:
    'someOutputFile in the respective *.zarr folder...'
  message: """--- Converting vcf into zarr format ---"""
  shell:
    """
    vcf2zarr.py {input.vcf} {output.    }
    """

