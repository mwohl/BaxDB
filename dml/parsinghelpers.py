import pandas as pd
import csv

# Create a list of chromosome names in the format "Chr#" for specified number of chromosomes
def generate_chromosome_list(numChromosomes):
  chrlist = []
  for count in range(1,numChromosomes+1):
    chrname = 'chr'+str(count)
    chrlist.append(chrname)
  return chrlist

def parse_lines_from_file(lineFile):
  linelist = []
  with open(lineFile) as f:
    rawlines = f.readlines()
    for linename in rawlines:
      linename = linename.rstrip()
      linelist.append(linename)
  return linelist

def convert_linelist_to_lineIDlist(conn, linelist, populationID):
  lineIDlist = []
  for linename in linelist:
    lineID = find.find_line(conn, linename, populationID)
    lineIDlist.append(lineID)
  return lineIDlist

def parse_variants_from_file(variantPosFile):
  with open(variantPosFile) as f:
    variantReader = csv.reader(f, delimiter='\t')
    variantlist = []
    for variant in variantReader:
      variantlist.append(variant[1])
  return variantlist

def parse_genotypes_from_file(genotypeFile):
  rawGenos = []
  with open(genotypeFile) as f:
    genoReader = csv.reader(f, delimiter='\t')
    for item in genoReader:
      # Remove first item, which is an index term
      item.pop(0)
      # Convert all genotype values to integers
      for i in range(len(item)):
        item[i] = int(item[i])
      rawGenos.append(item)
  return rawGenos

def parse_unique_runs_from_gwas_results_file(filepath):
  gwas_runs = []
  df = pd.read_csv(filepath)
  for index, row in df.iterrows():
    gwas_run = [row['trait'],row['nSNPs'],row['nLines']]
    if gwas_run not in gwas_runs:
      gwas_runs.append(gwas_run)
  return gwas_runs