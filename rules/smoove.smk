
rule smoove_call:
  input:
    bam = "bams/{sample}" + config['suffix'] + ".bam"
  output:
    out = "output/SMOOVE/indResult/{sample}-smoove.genotyped.vcf.gz"
  params:
    fasta = config["reference_fasta"],
    sample = "{sample}",
    outdir = "output/SMOOVE/indResult"
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
    expand("output/SMOOVE/indResult/{sample}-smoove.genotyped.vcf.gz", sample = samples)
  output:
    out = "output/SMOOVE/merged/smoove-merged.sites.vcf.gz"
  params:
    fasta = config['reference_fasta'],
    name = "smoove-merged",
    outdir = "output/SMOOVE/merged"
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

