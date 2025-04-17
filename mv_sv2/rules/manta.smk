
import pandas as pd

chrms = pd.read_csv(config['chr_list'], header = None)[0]
samples = pd.read_csv(config['manta_samples'], header = None)[0]

singularity: config['sif']

rule all:
  input:
    "../../svars/manta-merged/manta-merged.sites.vcf.gz"

rule manta_inv:
  input:
    vcf = "../../svars/manta/{sample}.manta.diploidSV.{ref}.vcf.gz"
  output:
    vcf = "../../svars/manta/{sample}.manta.inv.{ref}.vcf.gz"
  params:
    fasta = config['reference_fasta'],
    tmp = "output/MANTA/indResult/{sample}/{sample}.tmp.vcf.gz"
  threads: 1
  resources:
    time = 60,
    mem_mb = 20000
  shell:
    ''' 
      set +eu
      source activate manta
      set -e

      convertInversion.py \
        /opt/conda/bin/samtools \
        {params.fasta} \
        {input.vcf} |sed '/^\[/d' |bgzip -c > {output.vcf}

      tabix {output.vcf}
    '''

rule make_manta_list:
  input:
    expand("../../svars/manta/{sample}.manta.inv.{ref}.vcf.gz", sample = samples, ref = config['ref'])
  output:
    "../../lists/manta.vcf.list"
  shell:
    "ls {input} > {output}"

rule manta_svimmer:
  input:
    "../../lists/manta.vcf.list"
  output:
    vcf = "../../svars/manta-merged/separate/{chrm}.svimmer.vcf"
  params:
    chrm = "{chrm}" 
  threads: 6
  resources:
    time = 120,
    mem_mb = 60000
  shell:
    '''
      source activate svu

      ../../utils/svimmer/svimmer \
        --threads {threads} \
        {input} {params.chrm} > {output.vcf}
    '''

rule manta_join:
  input:
    expand("../../svars/manta-merged/separate/{chrm}.svimmer.vcf", chrm = chrms)
  output:
    vcf = "../../svars/manta-merged/manta-merged.sites.vcf.gz"
  threads: 2
  resources:
    time = 120,
    mem_mb = 20000
  shell:
    '''
      bcftools concat \
        -O z \
        -o {output.vcf} \
        {input}

      tabix {output.vcf}
    '''

