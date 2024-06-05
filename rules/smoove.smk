
rule smoove_call:
  input:
    bam = "bams/{sample}_pe_sorted_dedup.bam"
  output:
    out = "SMOOVE/indResult/{sample}-smoove.genotyped.vcf.gz"
  params:
    fasta = config["reference_fasta"],
    sample = "{sample}",
    outdir = "SMOOVE/indResult"
  threads: 4
  resources:
    time = 120,
    mem_mb = 60000
  shell:
    '''
      source activate smoove

      smoove call \
        --outdir {params.outdir} \
        --name {params.sample} \
        --fasta {params.fasta} \
        -p {threads} \
        --genotype {input.bam}
    '''

rule smoove_merge:
  input:
    expand("SMOOVE/indResult/{sample}-smoove.genotyped.vcf.gz", sample = samples)
  output:
    out = "SMOOVE/merged/smoove-merged.sites.vcf.gz"
  params:
    fasta = config['reference_fasta'],
    name = "smoove-merged",
    outdir = "SMOOVE/merged"
  threads: 20
  resources:
    time = 600,
    mem_mb = 60000
  shell:
    '''
      source activate smoove

      smoove merge \
        --name {params.name} \
        -f {params.fasta} \
        --outdir {params.outdir} \
        {input}
    '''

