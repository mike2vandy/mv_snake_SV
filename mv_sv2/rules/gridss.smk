
import pandas as pd

chrms = pd.read_csv(config['chr_list'], header = None)[0]
samples = pd.read_csv(config['gridss_samples'], header = None)[0]

singularity: config['sif']

rule all:
  input:
    "../../svars/gridss-merged/gridss-merged.sites.vcf.gz"

rule gridss_filter:
  input:
    vcf = "../../svars/gridss/{sample}.gridss.{ref}.vcf.gz"
  output:
    vcf = "../../svars/gridss/filter/{sample}.gridss.filter.{ref}.vcf"
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
    vcf = "../../svars/gridss/filter/{sample}.gridss.filter.{ref}.vcf"
  output:
    vcf = "../../svars/gridss/anno/{sample}.gridss.anno.{ref}.vcf.gz"
  params:
    tmp = "../../svars/gridss/anno/{sample}.gridss.anno.tmp.vcf"
  threads: 1
  resources:
    time = 120,
    mem_mb = 40000
  shell:
    '''
      set +eu
      source activate gridss
      set -e

      ../../utils/scripts/simple-event-annotation.R {input.vcf} {params.tmp}

      bcftools view \
        -e 'SVTYPE="CTX"' \
        -Oz \
        -o {output.vcf} \
        {params.tmp}

      tabix {output.vcf}
    '''

rule make_gridss_list:
  input:
    expand("../../svars/gridss/anno/{sample}.gridss.anno.{ref}.vcf.gz", sample = samples, ref = config['ref'])
  output:
    "../../lists/gridss.vcf.list"
  shell:
    "ls {input} > {output}"

rule gridss_svimmer:
  input:
    "../../lists/gridss.vcf.list"
  output:
    "../../svars/gridss-merged/separate/{chrm}.svimmer.vcf"
  params:
    chrm = "{chrm}"
  threads: 6
  resources:
    time = 120,
    mem_mb = 40000
  shell:
    '''
      ../../utils/svimmer/svimmer \
        --threads {threads} \
        {input} {params.chrm} > {output}
    '''

rule gridss_join:
  input:
    expand("../../svars/gridss-merged/separate/{chrm}.svimmer.vcf", chrm = chrms)
  output:
    "../../svars/gridss-merged/gridss-merged.sites.vcf.gz"
  threads: 1
  resources:
    time = 120,
    mem_mb = 4000
  shell:
    '''
      bcftools concat \
        -Oz \
        -o {output} \
        {input}
    '''

