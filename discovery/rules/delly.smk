
rule delly_merge:
  input:
    expand("/gpfs_common/share01/stern/mwvandew/make_SV_vcf/svars/delly/{sample}.delly.{ref}.vcf.gz", sample = delly_samples, ref = config['ref'])
  output:
    vcf = "output/delly/delly_merged/delly_merged.sites.vcf.gz"
  params:
    tmp = "output/delly/delly_merged/delly_merged.sites.bcf"
  threads: 1
  resources:
    time = 120,
    mem_mb = 60000
  shell:
    '''
    source activate delly

    delly merge {input} \
      -o {params.tmp}

    bcftools sort -Oz -o {output.vcf} {params.tmp}

    tabix {output.vcf}
    '''
