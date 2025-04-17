
import pandas as pd

samples = pd.read_csv(config['delly_samples'], header = None)[0]

singularity: config['sif']

rule all:
  input:
    "../../svars/delly-merged/delly-merged.sites.vcf.gz"

rule delly_merge:
  input:
    expand("../../svars/delly/{sample}.delly.{ref}.vcf.gz", sample = samples, ref = config['ref'])
  output:
    vcf = "../../svars/delly-merged/delly-merged.sites.vcf.gz"
  params:
    tmp = "../../svars/delly-merged/delly-merged.sites.bcf"
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
