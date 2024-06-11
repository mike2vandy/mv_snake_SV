
rule prep_gridss_ref:
  input:
    fasta = config['reference_fasta']
  output:
    bwt = config['reference_fasta'] + ".bwt"
  threads: 1
  resources:
    time = 240,
    mem_mb = 3600
  shell:
    '''
      set +eu
      source activate gridss
      set -e

      gridss \
        -s processreference \
        -r {input.fasta}
    '''

rule gridss_call:
  input:
    bam = "bams/{sample}" + config['suffix'] + ".bam",
    bwt = config['reference_fasta'] + ".bwt"
  output:
    bam = "output/GRIDSS/indResults/{sample}/{sample}.gridss.bam",
    vcf = "output/GRIDSS/indResults/{sample}/{sample}.gridss.vcf"
  params:
    fasta = config['reference_fasta'],
    workDir = "output/GRIDSS/indResults/{sample}"
  threads: 8
  resources:
    time = 240,
    mem_mb = 36000
  shell:
    '''
      set +eu
      source activate gridss
      set -e

      unset -f which 
    
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
    vcf = "output/GRIDSS/indResults/{sample}/{sample}.gridss.vcf"
  output:
    vcf = "output/GRIDSS/indResults/{sample}/{sample}.gridss.filter.vcf"
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
    vcf = "output/GRIDSS/indResults/{sample}/{sample}.gridss.filter.vcf"
  output:
    vcf = "output/GRIDSS/indResults/{sample}/{sample}.gridss.anno.vcf.gz"
  params:
    tmp = "output/GRIDSS/indResults/{sample}/{sample}.gridss.anno.tmp.vcf"
  threads: 1
  resources:
    time = 120,
    mem_mb = 40000
  shell:
    '''
      set +eu
      source activate gridss
      set -e

      utils/scripts/simple-event-annotation.R {input.vcf} {params.tmp}

      bcftools view \
        -e 'SVTYPE="CTX"' \
        -Oz \
        -o {output.vcf} \
        {params.tmp}

      tabix {output.vcf}
    '''

rule make_gridss_list:
  input:
    expand("output/GRIDSS/indResults/{sample}/{sample}.gridss.anno.vcf.gz", sample = samples)
  output:
    "output/lists/gridss.vcf.list"
  shell:
    "ls {input} > {output}"

rule gridss_svimmer:
  input:
    "output/lists/gridss.vcf.list"
  output:
    "output/GRIDSS/merged/separate/{chrm}.svimmer.vcf"
  params:
    chrm = "{chrm}"
  conda: "../env/pysam.yaml"
  threads: 8
  resources:
    time = 120,
    mem_mb = 40000
  shell:
    '''
      utils/svimmer/svimmer \
        --threads {threads} \
        {input} {params.chrm} > {output}
    '''

rule gridss_join:
  input:
    expand("output/GRIDSS/merged/separate/{chrm}.svimmer.vcf", chrm = chrms)
  output:
    "output/GRIDSS/merged/gridss-merged.sites.vcf.gz"
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

