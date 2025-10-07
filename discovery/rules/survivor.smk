
rule copyVCFs:
  input:
    "output/manta/manta_merged/manta_merged.sites.vcf.gz",
    "output/smoove/smoove_merged/smoove_merged.sites.vcf.gz",
    "output/delly/delly_merged/delly_merged.sites.vcf.gz",
    "output/gridss/gridss_merged/gridss_merged.sites.vcf.gz"
  output:
    "output/survivor/callers/manta_merged.sites.vcf",
    "output/survivor/callers/smoove_merged.sites.vcf",
    "output/survivor/callers/delly_merged.sites.vcf",
    "output/survivor/callers/gridss_merged.sites.vcf"
  threads: 1
  resources:
    time = 20,
    mem_mb = 20000
  shell:
    '''
      mkdir -p output/survivor/callers
      cp {input} output/survivor/callers
      gunzip output/survivor/callers/*
    '''

rule survivor:
  input:
    "output/survivor/callers/manta_merged.sites.vcf",
    "output/survivor/callers/smoove_merged.sites.vcf",
    "output/survivor/callers/delly_merged.sites.vcf",
    "output/survivor/callers/gridss_merged.sites.vcf"
  output:
    "output/survivor/merged/survivor.raw.vcf"
  threads: 1
  resources:
    time = 120,
    mem_mb = 60000
  shell:
    '''
      ls {input} > output/survivor/survivor.list

      source activate svu

      SURVIVOR merge output/survivor/survivor.list 100 1 1 0 0 50 {output}
    '''

rule filter_survivor:
  input:
    "output/survivor/merged/survivor.raw.vcf"
  output:
    vcf = "output/survivor/filter/survivor.filter.vcf" 
  params: 
    manta = "output/survivor/callers/manta_merged.sites.vcf",
  threads: 1
  resources:
    time = 60,
    mem_mb = 40000
  run:
    import vcf
    
    def format_line(record):
    
      info_fields = []
      for key, value in record.INFO.items():
        if isinstance(value, list):
          value_str = ','.join(map(str,value))
        else:
          value_str = str(value)
        info_fields.append('{}={}'.format(key, value_str))
    
      vcf_line = '\t'.join([
            record.CHROM,
            str(record.POS),
            record.ID if record.ID else '.',
            record.REF,
            ','.join(str(alt) for alt in record.ALT),
            str(record.QUAL) if record.QUAL is not None else '.',
            ','.join(record.FILTER) if record.FILTER else 'PASS',
            ';'.join(info_fields)])
    
      return vcf_line
    
    vcf_r = vcf.Reader(filename = str(input[0]))
    manta = params.manta
    
    with open(output.vcf, 'w') as out_vcf:
      for line in vcf_r._header_lines:
        print(line, file = out_vcf)
      print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO", file = out_vcf)
      
      for record in vcf_r:
        info = record.INFO
        if info['SVTYPE'] != 'TRA':
          if info['SUPP'] == '1' and info['SVTYPE'] == 'INS':
            manta_geno = record.genotype(manta)
            if manta_geno['TY'] == 'INS':
              print(format_line(record), file = out_vcf)
          elif int(info['SUPP']) >= 2:
             print(format_line(record), file = out_vcf) 

rule bgzip_vcf:
  input:
    "output/survivor/filter/survivor.filter.vcf"
  output:
    vcf = "output/survivor/filter/survivor.filter.vcf.gz",
    tbi =  "output/survivor/filter/survivor.filter.vcf.gz.tbi"
  threads: 1
  resources:
    time = 60,
    mem_mb = 40000
  shell:
    '''
      bcftools sort -Oz -o {output.vcf} {input}
      tabix {output.vcf} 
    '''
