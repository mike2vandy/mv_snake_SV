
rule manta_call:
  input:
    bam = "bams/{sample}_pe_sorted_dedup.bam"
  output:
    out = "MANTA/indResult/{sample}/results/variants/diploidSV.vcf.gz"
  params:
    fasta = "ref/rheMac10.fa",
    workDir = "MANTA/indResult/{sample}"
  threads: 24
  resources:
    time = 600,
    mem_mb = 60000
  shell:
    '''
      set +eu
      source activate manta
      set -e

      configManta.py \
        --bam {input.bam} \
        --referenceFasta {params.fasta} \
        --runDir {params.workDir}

      cd {params.workDir}

      ./runWorkflow.py \
        --quiet \
        -m local \
        -j {threads}
    '''

rule manta_inv:
  input:
    vcf = "MANTA/indResult/{sample}/results/variants/diploidSV.vcf.gz"
  output:
    vcf = "MANTA/indResult/{sample}/{sample}.manta.vcf.gz"
  params:
    fasta = "ref/rheMac10.fa",
    tmp = "MANTA/indResult/{sample}/{sample}.tmp.vcf.gz"
  threads: 1
  resources:
    time = 60,
    mem_mb = 20000
  shell:
    '''
      bcftools view \
        -f 'PASS' \
        -Oz \
        -o {params.tmp} \
        {input.vcf}

      set +eu
      source activate manta
      set -e

      convertInversion.py \
        /opt/conda/bin/samtools \
        {params.fasta} \
        {params.tmp} |sed '/^\[/d' |bgzip -c > {output.vcf}

      tabix {output.vcf}
    '''

rule make_manta_list:
  input:
    expand("MANTA/indResult/{sample}/{sample}.manta.vcf.gz", sample = samples)
  output:
    "manta.vcf.list"
  shell:
    "ls {input} > {output}"

rule manta_svimmer:
  input:
    "manta.vcf.list"
  output:
    vcf = "MANTA/groupCall/separate/{chrm}.svimmer.vcf"
  params:
    chrm = "{chrm}"
  conda: "env/pysam.yaml"
  threads: 6
  resources:
    time = 600,
    mem_mb = 60000
  shell:
    '''
      ./svimmer/svimmer \
        --threads {threads} \
        {input} {params.chrm} > {output.vcf}
    '''

rule manta_join:
  input:
    expand("MANTA/groupCall/separate/{chrm}.svimmer.vcf", chrm = chrms)
  output:
    vcf = "MANTA/groupCall/manta-merged.sites.vcf.gz"
  threads: 2
  resources:
    time = 120,
    mem_mb = 60000
  shell:
    '''
      bcftools concat \
        -O z \
        -o {output.vcf} \
        {input}

      tabix {output.vcf}
    '''

