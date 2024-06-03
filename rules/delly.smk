
rule delly_call:
  input:
    bam = "bams/{sample}_pe_sorted_dedup.bam"
  output:
    vcf = "DELLY/indResults/{sample}/filter/{sample}_{sv_type}.vcf.gz"
  params:
    sv = '{sv_type}',
    fasta = "ref/rheMac10.fa",
    tmp = "DELLY/indResults/{sample}/raw/{sample}_{sv_type}.bcf"
  threads: 1
  resources:
    time = 600,
    mem_mb = 60000
  shell:
    '''
      source activate delly

      mkdir -p DELLY/indResults/{wildcards.sample}/raw

      delly call \
        -t {params.sv} \
        -g {params.fasta} \
        -o {params.tmp} \
        {input.bam}

      bcftools view \
        -O z \
        -o {output.vcf} \
        -f 'PASS' \
        {params.tmp}

      tabix {output.vcf}
    '''

rule delly_concat:
  input:
    "DELLY/indResults/{sample}/filter/{sample}_DEL.vcf.gz",
    "DELLY/indResults/{sample}/filter/{sample}_DUP.vcf.gz",
    "DELLY/indResults/{sample}/filter/{sample}_INV.vcf.gz",
    "DELLY/indResults/{sample}/filter/{sample}_INS.vcf.gz"
  output:
    vcf = "DELLY/indResults/{sample}/{sample}.delly.vcf.gz"
  threads: 1
  resources:
    time = 120,
    mem_mb = 60000
  shell:
    '''
      bcftools concat \
        -a -Ov \
        -o {output.vcf} \
        {input}

      tabix {output.vcf}
    '''

rule delly_merge:
  input:
    expand("DELLY/indResults/{sample}/{sample}.delly.vcf.gz", sample = samples)
  output:
    vcf = "DELLY/merged/delly-merged.sites.vcf.gz"
  params:
    tmp = "DELLY/merged/delly-merged.sites.bcf"
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



