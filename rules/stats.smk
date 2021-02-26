rule vcf2zarr:
  input:
    vcf = 'vars/ta{vartype}.vcf'
  output:
    'vars/ta{vartype}/.zgroup'
  message: """--- Converting vcf into zarr format ---"""
  shell:
    """
    echo 'hello'
##vcf2zarr.py {input.vcf} {output.    }
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


