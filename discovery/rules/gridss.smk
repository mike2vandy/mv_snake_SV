
rule gridss_filter:
  input:
    vcf = "/gpfs_common/share01/stern/mwvandew/make_SV_vcf/svars/gridss/{sample}.gridss.{ref}.vcf.gz"
  output:
    vcf = "output/gridss/filter/{sample}.gridss.filter.{ref}.vcf"
  threads: 1
  resources:
    time = 120,
    mem_mb = 40000
  shell:
    '''
      bcftools view -f "PASS" {input.vcf} > {output.vcf}

      sed -i '100i\##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="SV_length">' {output.vcf}
    '''

rule gridss_anno:
  input:
    vcf = "output/gridss/filter/{sample}.gridss.filter.{ref}.vcf"
  output:
    vcf = "output/gridss/anno/{sample}.gridss.anno.{ref}.vcf.gz"
  params:
    tmp = "output/gridss/anno/{sample}.gridss.anno.tmp.vcf"
  threads: 1
  resources:
    time = 120,
    mem_mb = 40000
  shell:
    '''
      set +eu
      source activate gridss
      set -e

      ./utils/scripts/simple-event-annotation.R {input.vcf} {params.tmp}

      bcftools view \
        -e 'SVTYPE="CTX"' \
        -Oz \
        -o {output.vcf} \
        {params.tmp}

      tabix {output.vcf}
    '''

rule make_gridss_list:
  input:
    expand("output/gridss/anno/{sample}.gridss.anno.{ref}.vcf.gz", sample = gridss_samples, ref = config['ref'])
  output:
    "output/vcf_lists/gridss.vcf.list"
  shell:
    "ls {input} > {output}"

rule gridss_svimmer:
  input:
    "output/vcf_lists/gridss.vcf.list"
  output:
    "output/gridss/gridss_merged/separate/{chrm}.svimmer.vcf"
  params:
    chrm = "{chrm}"
  threads: 6
  resources:
    time = 240,
    mem_mb = 40000
  shell:
    '''
      source activate svu

      /gpfs_common/share01/stern/mwvandew/make_SV_vcf/svimmer/svimmer \
        --threads {threads} \
        {input} {params.chrm} > {output}
    '''

rule gridss_join:
  input:
    expand("output/gridss/gridss_merged/separate/{chrm}.svimmer.vcf", chrm = chrms)
  output:
    "output/gridss/gridss_merged/gridss_merged.sites.vcf.gz"
  threads: 1
  resources:
    time = 120,
    mem_mb = 20000
  shell:
    '''
      bcftools concat \
        -Oz \
        -o {output} \
        {input}
    '''

