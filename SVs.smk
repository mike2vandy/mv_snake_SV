import pandas as pd

#read sample file list 
tmpPd = pd.read_table(config['sample_list'], header = None)
samples = list(tmpPd[0])

#read chromosome file list
tmpCh = pd.read_table(config['chr_list'], header = None)
chrms = list(tmpCh[0])

#svs for delly, no bend
delSVs = ['DEL', 'INS', 'DUP', 'INV', 'BND']

#path to singularity container 
singularity: config['sif']

rule all:
  input:
    "SMOOVE/merged/smoove-merged.sites.vcf.gz", 
    expand("GRAPHTYPER/final/{species}.{reference}.SVs.vcf.gz", 
      species = config['species'], 
      reference = config['reference']), 
    #expand("GRAPHTYPER/{sample}/{sample}.genotype.vcf.gz", sample = samples)
    "SURVIVOR/filter/survivor.filter.vcf.gz",
    #"SURVIVOR/main/survivor.raw.vcf"
    "MANTA/merged/manta-merged.sites.vcf.gz",
    "DELLY/merged/delly-merged.sites.vcf.gz",
    "GRIDSS/merged/gridss-merged.sites.vcf.gz"
    #expand("GRIDSS/indResults/{sample}/{sample}.gridss.vcf", sample = samples)
    #expand("DELLY/indResults/{sample}/filter/{sample}_{sv_type}.vcf.gz", sample = samples, sv_type = delSVs)

include: "rules/smoove.smk"
include: "rules/manta.smk"
include: "rules/delly.smk"
include: "rules/gridss.smk"
include: "rules/survivor.smk"
include: "rules/graph-typer.smk"


