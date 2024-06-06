
rule copyVCFs:
  input:
    "output/MANTA/merged/manta-merged.sites.vcf.gz",
    "output/SMOOVE/merged/smoove-merged.sites.vcf.gz",
    "output/DELLY/merged/delly-merged.sites.vcf.gz",
    "output/GRIDSS/merged/gridss-merged.sites.vcf.gz"
  output:
    "output/SURVIVOR/callers/manta-merged.sites.vcf",
    "output/SURVIVOR/callers/smoove-merged.sites.vcf",
    "output/SURVIVOR/callers/delly-merged.sites.vcf",
    "output/SURVIVOR/callers/gridss-merged.sites.vcf"
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
    "output/SURVIVOR/callers/manta-merged.sites.vcf",
    "output/SURVIVOR/callers/smoove-merged.sites.vcf",
    "output/SURVIVOR/callers/delly-merged.sites.vcf",
    "output/SURVIVOR/callers/gridss-merged.sites.vcf"
  output:
    "output/SURVIVOR/main/survivor.raw.vcf"
  threads: 1
  resources:
    time = 120,
    mem_mb = 60000
  shell:
    '''
      ls {input} > output/lists/survivor.list

      source activate svu

      SURVIVOR merge output/lists/survivor.list 100 1 1 0 0 50 {output}
    '''

rule filter_survivor:
  input:
    "output/SURVIVOR/main/survivor.raw.vcf"
  output:
    "output/SURVIVOR/filter/survivor.filter.vcf.gz"
  params:
    tmp = "output/SURVIVOR/main/survivor.filter.vcf"
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

