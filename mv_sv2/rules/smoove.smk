
import pandas as pd

singularity: config['sif']

samples = pd.read_csv(config['manta_samples'], header = None)[0]

rule all:
  input:
    "../../svars/smoove_merged/smoove-merged.sites.vcf.gz"

rule smoove_merge:
  input:
    expand("../../svars/smoove/{sample}.smoove.{ref}.vcf.gz", sample = samples, ref = config['ref'])
  output:
    out = "../../svars/smoove_merged/smoove-merged.sites.vcf.gz"
  params:
    fasta = config['reference_fasta'],
    name = "smoove-merged",
    outdir = "../../svars/smoove_merged"
  threads: 12
  resources:
    time = 120,
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

