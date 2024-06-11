#! /usr/bin/env python

import vcf, sys

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

vcfR = vcf.Reader(filename = sys.argv[1])

manta = "output/SURVIVOR/callers/manta-merged.sites.vcf"

for line in vcfR._header_lines:
  print(line)
print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO")

for record in vcfR:
  info = record.INFO
  if info['SVTYPE'] != 'TRA':
    if info['SUPP'] == '1' and info['SVTYPE'] == 'INS':
      manta_geno = record.genotype(manta)
      if manta_geno['TY'] == 'INS':
        print(format_line(record))
    elif int(info['SUPP']) >= 2:
       print(format_line(record))

