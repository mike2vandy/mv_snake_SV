
rule gridss_call:
  input:
    bam = "bams/{sample}_pe_sorted_dedup.bam"
  output:
    bam = "GRIDSS/indResults/{sample}/{sample}.gridss.bam",
    vcf = "GRIDSS/indResults/{sample}/{sample}.gridss.vcf"
  params:
    fasta = "ref/rheMac10.fa",
    workDir = "GRIDSS/indResults/{sample}"
  threads: 8
  resources:
    time = 240,
    mem_mb = 36000
  shell:
    '''
      set +eu
      source activate gridss
      set -e

      gridss \
        -t {threads} \
        -r {params.fasta} \
        -o {output.vcf} \
        -a {output.bam} \
        --jvmheap 32g \
        -w {params.workDir} \
        {input.bam}
    '''

rule gridss_filter:
  input:
    vcf = "GRIDSS/indResults/{sample}/{sample}.gridss.vcf"
  output:
    vcf = "GRIDSS/indResults/{sample}/{sample}.gridss.filter.vcf"
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
    vcf = "GRIDSS/indResults/{sample}/{sample}.gridss.filter.vcf"
  output:
    vcf = "GRIDSS/indResults/{sample}/{sample}.gridss.anno.vcf.gz"
  params:
    tmp = "GRIDSS/indResults/{sample}/{sample}.gridss.anno.tmp.vcf"
  threads: 1
  resources:
    time = 120,
    mem_mb = 40000
  shell:
    '''
      set +eu
      source activate gridss
      set -e

      ./software/simple-event-annotation.R {input.vcf} {params.tmp}

      bcftools view \
        -e 'SVTYPE="CTX"' \
        -Oz \
        -o {output.vcf} \
        {params.tmp}

      tabix {output.vcf}
    '''

rule make_gridss_list:
  input:
    expand("GRIDSS/indResults/{sample}/{sample}.gridss.anno.vcf.gz", sample = samples)
  output:
    "gridss.vcf.list"
  shell:
    "ls {input} > {output}"

rule gridss_svimmer:
  input:
    "gridss.vcf.list"
  output:
    "GRIDSS/groupCall/separate/{chrm}.svimmer.vcf"
  params:
    chrm = "{chrm}"
  conda: "env/pysam.yaml"
  threads: 8
  resources:
    time = 120,
    mem_mb = 40000
  shell:
    '''
      ~/software/svimmer/svimmer \
        --threads {threads} \
        {input} {params.chrm} > {output}
    '''

rule gridss_join:
  input:
    expand("GRIDSS/groupCall/separate/{chrm}.svimmer.vcf", chrm = chrms)
  output:
    "GRIDSS/groupCall/gridss-merged.sites.vcf.gz"
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

