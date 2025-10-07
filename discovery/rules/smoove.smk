
 
rule smoove_merge:
  input:
    expand("/gpfs_common/share01/stern/mwvandew/make_SV_vcf/svars/smoove/{sample}.smoove.{ref}.vcf.gz", sample = smoove_samples, ref = config['ref'])
  output:
    out = "output/smoove/smoove_merged/smoove_merged.sites.vcf.gz"
  params:
    fasta = config['reference_fasta'],
    name = "smoove_merged",
    outdir = "output/smoove/smoove_merged"
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

