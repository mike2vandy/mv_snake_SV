
rule graphtyper:
  input:
    cram = f"/gpfs_common/share01/stern/mwvandew/make_SV_vcf/crams/{{sample}}.{config['suffix']}.cram",
    vcf = "output/survivor/filter/survivor.filter.vcf.gz"
  output:
    vcf = "output/graphtyper/individuals/{sample}/{sample}.genotype.vcf.gz"
  params:
    out_dir = "output/graphtyper/individuals/{sample}/chrms",
    fasta = config['reference_fasta'],
    region = config['chr_list'], 
    cache = config['cache']
  threads: 10
  resources:
    time = 360,
    mem_mb = 40000
  shell:
    '''
      source activate svu
      
      export REF_PATH={params.cache}/%2s/%2s/%s
      export REF_CACHE={params.cache}/%2s/%2s/%s

      graphtyper genotype_sv \
        {params.fasta} \
        {input.vcf} \
        --force_no_filter_zero_qual \
        --threads {threads} \
        --sam {input.cram} \
        --region_file {params.region} \
        --output {params.out_dir}

      set +e

      cat {params.region} \
        |while read chrom; do \
           if [[ ! -d {params.out_dir}/${{chrom}} ]]; then continue; fi; \
           find {params.out_dir}/${{chrom}} -name "*.vcf.gz" | sort ; \
         done > {params.out_dir}/vcf.list

      bcftools concat -n -f {params.out_dir}/vcf.list \
        -Oz -o {output.vcf}

      tabix {output.vcf}
    '''

rule graph_merge:
  input:
    expand("output/graphtyper/individuals/{sample}/{sample}.genotype.vcf.gz", sample = cram_samples)
  output:
    vcf = f"output/graphtyper/main_raw/{config['species']}.square.vcf.gz"
  threads: 1
  resources:
    time = 240,
    mem_mb = 60000
  shell:
    '''
      source activate svu 

      graphtyper vcf_merge \
        --sv \
        {input} |bgzip -c > {output.vcf}

      tabix {output.vcf}
    '''

rule graph_filter:
  input:
    f"output/graphtyper/main_raw/{config['species']}.square.vcf.gz"
  output:
    f"output/graphtyper/main_filter/{config['species']}.{config['ref']}.SVs.vcf.gz"
  threads: 1
  resources:
    time = 60,
    mem_mb = 40000
  shell:
    '''
      bcftools view \
        -f "PASS" \
        -i "SVMODEL='AGGREGATED'" \
        -Oz \
        -o {output} \
        {input}

      tabix {output}
    '''

rule graph_vep:
  input:
    vcf = f"output/graphtyper/main_filter/{config['species']}.{config['ref']}.SVs.vcf.gz"
  output:
    out_tmp = temp(f"output/graphtyper/main_filter/{config['species']}.{config['ref']}.SVs.vep.vcf"),
    final = f"output/graphtyper/main_filter/{config['species']}.{config['ref']}.SVs.vep.vcf.gz"
  params:
    fasta = config['reference_fasta'],
    gtf = config['gtf']
  threads: 4
  resources:
    time = 120,
    mem_mb = 60000
  shell:
    '''
      set -e
      source activate ensembl-vep

      vep \
        -i {input.vcf} \
        -o {output.out_tmp} \
        -gtf {params.gtf} \
        --fasta {params.fasta} \
        --format vcf \
        --vcf \
        --everything \
        --dont_skip \
        --fork {threads}

      bgzip -c {output.out_tmp} > {output.final}

      tabix {output.final}
    '''

