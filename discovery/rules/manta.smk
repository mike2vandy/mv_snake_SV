
import pandas as pd

rule manta_inv:
  input:
    vcf = "/gpfs_common/share01/stern/mwvandew/make_SV_vcf/svars/manta/{sample}.manta.diploidSV.{ref}.vcf.gz"
  output:
    vcf = "output/manta/inv/{sample}.manta.inv.{ref}.vcf.gz"
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
    expand("output/manta/inv/{sample}.manta.inv.{ref}.vcf.gz", sample = manta_samples, ref = config['ref'])
  output:
    "output/vcf_lists/manta.vcf.list"
  shell:
    "ls {input} > {output}"

rule manta_svimmer:
  input:
    "output/vcf_lists/manta.vcf.list"
  output:
    vcf = "output/manta/manta_merged/separate/{chrm}.svimmer.vcf"
  params:
    chrm = "{chrm}" 
  threads: 6
  resources:
    time = 120,
    mem_mb = 60000
  shell:
    '''
      source activate svu

      /gpfs_common/share01/stern/mwvandew/make_SV_vcf/svimmer/svimmer \
        --threads {threads} \
        {input} {params.chrm} > {output.vcf}
    '''

rule manta_join:
  input:
    expand("output/manta/manta_merged/separate/{chrm}.svimmer.vcf", chrm = chrms)
  output:
    vcf = "output/manta/manta_merged/manta_merged.sites.vcf.gz"
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

