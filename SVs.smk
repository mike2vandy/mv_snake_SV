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
    expand("output/GRAPHTYPER/final/{species}.{reference}.SVs.vep.vcf.gz",
      species = config['species'],
      reference = config['reference']) 

include: "rules/smoove.smk"
include: "rules/manta.smk"
include: "rules/delly.smk"
include: "rules/gridss.smk"
include: "rules/survivor.smk"
include: "rules/graph-typer.smk"


