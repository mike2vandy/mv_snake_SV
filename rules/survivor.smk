
rule copyVCFs:
  input:
    "MANTA/merged/manta-merged.sites.vcf.gz",
    "SMOOVE/merged/smoove-merged.sites.vcf.gz",
    "DELLY/merged/delly-merged.sites.vcf.gz",
    "GRIDSS/merged/gridss-merged.sites.vcf.gz"
  output:
    "SURVIVOR/callers/manta-merged.sites.vcf",
    "SURVIVOR/callers/smoove-merged.sites.vcf",
    "SURVIVOR/callers/delly-merged.sites.vcf",
    "SURVIVOR/callers/gridss-merged.sites.vcf"
  threads: 1
  resources:
    time = 20,
    mem_mb = 20000
  shell:
    '''
      cp {input} SURVIVOR/callers
      gunzip SURVIVOR/callers/*
    '''

rule survivor:
  input:
    "SURVIVOR/callers/manta-merged.sites.vcf",
    "SURVIVOR/callers/smoove-merged.sites.vcf",
    "SURVIVOR/callers/delly-merged.sites.vcf",
    "SURVIVOR/callers/gridss-merged.sites.vcf"
  output:
    "SURVIVOR/main/survivor.raw.vcf"
  threads: 1
  resources:
    time = 120,
    mem_mb = 60000
  shell:
    '''
      ls {input} > survivor.list

      source activate svu

      SURVIVOR merge survivor.list 100 1 1 0 0 50 {output}
    '''

rule filter_survivor:
  input:
    "SURVIVOR/main/survivor.raw.vcf"
  output:
    "SURVIVOR/filter/survivor.filter.vcf.gz"
  params:
    tmp = "SURVIVOR/main/survivor.filter.vcf"
  conda: "../env/gter.yaml"
  threads: 1
  resources:
    time = 60,
    mem_mb = 40000
  shell:
    '''
      utils/scripts/filterSurv.py {input} > {params.tmp}

      bcftools sort -Oz -o {output} {params.tmp}
      
      sleep 60

      /opt/conda/bin/tabix {output}
    '''

