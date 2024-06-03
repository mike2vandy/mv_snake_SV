import pandas as pd

monkeys = pd.read_table("samples.txt")
samples = list(monkeys['monkId'])
tmp = pd.read_table("ref/rhe.mac.chrm.list")
chrms = list(tmp['chrm'])

#svs for delly, no bend
delSVs = ['DEL', 'INS', 'DUP', 'INV']

#path to singularity container 
singularity: "/home/fried255/mvandewe/universalDat/UU_Cfam_GSD_ROSY/wags.sif"

rule all:
  input:
    "GRAPHTYPER/final/rheMac.square.filter.vcf.gz" 
    #expand("GRAPHTYPER/{sample}/{sample}.genotype.vcf.gz", sample = samples)
    #"SURVIVOR/filter/survivor.filter.vcf.gz"
    #"SURVIVOR/main/survivor.raw.vcf"
    #"MANTA/groupCall/manta-merged.sites.vcf.gz",
    #"SMOOVE/merged/smoove-merged.sites.vcf.gz",
    #"DELLY/merged/delly-merged.sites.vcf.gz",
    #"GRIDSS/groupCall/gridss-merged.sites.vcf.gz"
    #expand("GRIDSS/indResults/{sample}/{sample}.gridss.vcf", sample = samples)
    #expand("DELLY/indResults/{sample}/filter/{sample}_{sv_type}.vcf.gz", sample = samples, sv_type = delSVs)
 
