

rule graphTyper:
  input:
    bam = "bams/{sample}_pe_sorted_dedup.bam",
    vcf = "SURVIVOR/filter/survivor.filter.vcf.gz"
  output:
    vcf = "GRAPHTYPER/{sample}/{sample}.genotype.vcf.gz"
  params:
    outDir = "GRAPHTYPER/{sample}",
    fasta = config['reference_fasta'],
    region = config['chr_list']
  conda: "../env/gter.yaml"
  threads: 10
  resources:
    time = 240,
    mem_mb = 40000
  shell:
    '''
      graphtyper genotype_sv \
        {params.fasta} \
        {input.vcf} \
        --force_no_filter_zero_qual \
        --threads {threads} \
        --sam {input.bam} \
        --region_file {params.region} \
        --output {params.outDir}

      set +e

      cat {params.region} \
        |while read chrom; do if [[ ! -d {params.outDir}/${{chrom}} ]]; then continue; fi; find {params.outDir}/${{chrom}} -name "*.vcf.gz" | sort ; done > {params.outDir}/vcf.list

      bcftools concat -n -f {params.outDir}/vcf.list \
        -Oz -o {output.vcf}

      tabix {output.vcf}

    '''

rule graphMerge:
  input:
    expand("GRAPHTYPER/{sample}/{sample}.genotype.vcf.gz", sample = samples)
  output:
    vcf = expand("GRAPHTYPER/final/{species}.square.vcf.gz", species = config['species'])
  conda: "../env/gter.yaml"
  threads: 1
  resources:
    time = 60,
    mem_mb = 60000
  shell:
    '''
      graphtyper vcf_merge \
        --sv \
        {input} |bgzip -c > {output.vcf}

      tabix {output.vcf}
    '''

rule graphFilter:
  input:
    expand("GRAPHTYPER/final/{species}.square.vcf.gz", species = config['species'])
  output:
    expand("GRAPHTYPER/final/{species}.{reference}.SVs.vcf.gz", species = config['species'], reference = config['reference'])
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

