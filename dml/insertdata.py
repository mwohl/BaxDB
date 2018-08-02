from configparser import ConfigParser
import pandas as pd
import psycopg2
import csv

def config(filename='database.ini', section='postgresql'):
  # create a parser
  parser = ConfigParser()
  # read config file
  parser.read(filename)

  # get section, default to postgresql
  db = {}
  if parser.has_section(section):
    params = parser.items(section)
    for param in params:
      db[param[0]] = param[1]
  else:
    raise Exception('Section {0} not found in the {1} file'.format(section, filename))

  return db

def connect():
  conn = None
  try:
    # read connection parameters
    params = config()

    # connect to the PostgreSQL server
    print('Connection to the PostgreSQL database...')
    conn = psycopg2.connect(**params)

    # create a cursor
    cur = conn.cursor()

    print('PostgreSQL database version:')
    cur.execute('SELECT version()')

    db_version = cur.fetchone()
    print(db_version)

    # Close the connection
    cur.close()
  except (Exception, psycopg2.DatabaseError) as error:
    print(error)
  
  return conn

def insert_species(conn, species):
  cur = conn.cursor()
  SQL = """INSERT INTO species (shortname, binomial, subspecies, variety)
        VALUES (%s, %s, %s, %s)
        ON CONFLICT DO NOTHING
        RETURNING species_id;"""
  args_tuple = (species.n, species.b, species.s, species.v)
  cur.execute(SQL, args_tuple)
  row = cur.fetchone()
  if row is not None:
    newID = row[0]
    conn.commit()
    cur.close()
    return newID
  else:
    return None

def find_species(conn, speciesShortname):
  cur = conn.cursor()
  cur.execute("SELECT species_id FROM species WHERE shortname = '%s';" %  speciesShortname)
  row = cur.fetchone()
  if row is not None:
    speciesID = row[0]  
    cur.close()
    return speciesID
  else:
    return None

def insert_population(conn, population): 
  cur = conn.cursor()
  SQL = """INSERT INTO population (population_name, population_species)
        VALUES (%s, %s)
        ON CONFLICT DO NOTHING
        RETURNING population_id;"""
  args_tuple = (population.n, population.s)
  cur.execute(SQL, args_tuple)
  row = cur.fetchone()
  if row is not None:
     newID = row[0]
     conn.commit()
     cur.close()
     return newID
  else:
    return None  

def find_population(conn, populationName):
  cur = conn.cursor()
  cur.execute("SELECT population_id FROM population WHERE population_name = '%s';" % populationName)
  row = cur.fetchone()
  if row is not None:
    populationID = row[0]
    cur.close()
    return populationID
  else:
    return None

def generate_chromosome_list(numChromosomes):
  chrlist = []
  for count in range(1,numChromosomes+1):
    chrname = 'chr'+str(count)
    chrlist.append(chrname)
  return chrlist

def insert_chromosome(conn, chromosome):
  cur = conn.cursor()
  SQL = """INSERT INTO chromosome (chromosome_name, chromosome_species)
        VALUES (%s, %s)
        ON CONFLICT DO NOTHING
        RETURNING chromosome_id;"""
  args_tuple = (chromosome.n, chromosome.s)
  cur.execute(SQL, args_tuple)
  row = cur.fetchone()
  if row is not None:
    newID = row[0]
    conn.commit()
    cur.close()
    return newID
  else:
    return None

def insert_all_chromosomes_for_species(conn, numChromosomes, speciesID):
  chrlist = generate_chromosome_list(numChromosomes)
  insertedChromosomeIDs = []
  for chrname in chrlist:
    chrobj = chromosome(chrname, speciesID)
    insertedChromosomeID = insert_chromosome(conn, chrobj)
    insertedChromosomeIDs.append(insertedChromosomeID)
  return insertedChromosomeIDs

def find_chromosome(conn, chromosome_name, chromosome_species):
  cur = conn.cursor()
  # not sure if next line is correct...
  cur.execute("SELECT chromosome_id FROM chromosome WHERE chromosome_name = '%s' AND chromosome_species = '%s';" % (chromosome_name, chromosome_species))
  row = cur.fetchone()
  if row is not None:
    chromosomeID = row[0]
    cur.close()
    return chromosomeID
  else:
    return None

def parse_lines_from_file(lineFile):
  linelist = []
  with open(lineFile) as f:
    rawlines = f.readlines()
    for linename in rawlines:
      linename = linename.rstrip()
      linelist.append(linename)
  return linelist

def insert_line(conn, line):
  cur = conn.cursor()
  SQL = """INSERT INTO line (line_name, line_population)
        VALUES (%s, %s)
        ON CONFLICT DO NOTHING
        RETURNING line_id;"""
  args_tuple = (line.n, line.p)
  cur.execute(SQL, args_tuple)
  row = cur.fetchone()
  if row is not None:
    newID = row[0]
    conn.commit()
    cur.close()
    return newID
  else:
    return None

def insert_lines_from_file(conn, lineFile, populationID):
  linelist = parse_lines_from_file(lineFile)
  insertedLineIDs = []
  for linename in linelist:
    lineobj = line(linename, populationID)
    insertedLineID = insert_line(conn, lineobj)
    insertedLineIDs.append(insertedLineID)
  return insertedLineIDs

def find_line(conn, line_name, line_population):
  cur = conn.cursor()
  cur.execute("SELECT line_id FROM line WHERE line_name = '%s' AND line_population = '%s';" % (line_name, line_population))
  row = cur.fetchone()
  if row is not None:
    lineID = row[0]
    cur.close()
    return lineID
  else:
    return None

def convert_linelist_to_lineIDlist(conn, linelist, populationID):
  lineIDlist = []
  for linename in linelist:
    lineID = find_line(conn, linename, populationID)
    lineIDlist.append(lineID)
  return lineIDlist

def parse_variants_from_file(variantPosFile):
  with open(variantPosFile) as f:
    variantReader = csv.reader(f, delimiter='\t')
    variantlist = []
    for variant in variantReader:
      variantlist.append(variant[1])
  return variantlist

def insert_variant(conn, variant):
  cur = conn.cursor()
  SQL = """INSERT INTO variant(variant_species, variant_chromosome, variant_pos)
        VALUES (%s,%s,%s)
        ON CONFLICT DO NOTHING
        RETURNING variant_id;"""
  args_tuple = (variant.s, variant.c, variant.p)
  cur.execute(SQL, args_tuple)
  #newID = cur.fetchone()[0]
  row = cur.fetchone()
  if row is not None:
    newID = row[0]
    conn.commit()
    cur.close()
    return newID
  else:
    return None

def insert_variants_from_file(conn, variantPosFile, speciesID, chromosomeID):
  variantlist = parse_variants_from_file(variantPosFile)
  print('num variants:')
  print(len(variantlist))
  insertedVariantIDs = []
  for variantpos in variantlist:
    variantobj = variant(speciesID, chromosomeID, variantpos)
    insertedVariantID = insert_variant(conn, variantobj)
    insertedVariantIDs.append(insertedVariantID)
  return insertedVariantIDs

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

def insert_genotypes_from_file(conn, genotypeFile, lineFile, chromosomeID, populationID):
  genotypes = parse_genotypes_from_file(genotypeFile)
  linelist = parse_lines_from_file(lineFile)
  lineIDlist = convert_linelist_to_lineIDlist(conn, linelist, populationID)
  zipped = zip(lineIDlist, genotypes)
  ziplist = list(zipped)
  insertedGenotypeIDs = []
  for zippedpair in ziplist:
    genotypeObj = genotype(zippedpair[0], chromosomeID, zippedpair[1])
    insertedGenotypeID = insert_genotype(conn, genotypeObj)
    insertedGenotypeIDs.append(insertedGenotypeID)
  return insertedGenotypeIDs

def insert_genotype(conn, genotype):
  cur = conn.cursor()
  SQL = """INSERT INTO genotype(genotype_line, genotype_chromosome, genotype)
        VALUES (%s,%s,%s)
        ON CONFLICT DO NOTHING
        RETURNING genotype_id;"""
  args_tuple = (genotype.l, genotype.c, genotype.g)
  cur.execute(SQL, args_tuple)
  newID = cur.fetchone()[0]
  conn.commit()
  cur.close()
  return newID

def insert_trait(conn, trait):
  cur = conn.cursor()
  SQL = """INSERT INTO trait(trait_name)
        VALUES (%s)
        ON CONFLICT DO NOTHING
        RETURNING trait_id;"""
  arg = (trait.n,)
  cur.execute(SQL, arg)
  row = cur.fetchone()
  if row is not None:
    newID = row[0]
    conn.commit()
    cur.close()
    return newID
  else:
    return None

def find_trait(conn, traitName):
  cur = conn.cursor()
  cur.execute("SELECT trait_id FROM trait WHERE trait_name = '%s';" % traitName)
  row = cur.fetchone()
  if row is not None:
    traitID = row[0]
    cur.close()
    return traitID
  else:
    return None

def insert_phenotypes_from_file(conn, phenotypeFile, populationID):
  phenotypeRawData = pd.read_csv(phenotypeFile, index_col=0)
  insertedPhenoIDs = []
  for key, value in phenotypeRawData.iteritems():
    print("***********KEY**************:")
    print key
    traitID = find_trait(conn, key)
    for index, traitval in value.iteritems():
      print("index:")
      print(index)
      lineID = find_line(conn, index, maize282popID)
      if lineID is None:
        newline = line(index, maize282popID)
        lineID = insert_line(conn, newline)
      print("trait value:")
      print(traitval)
      pheno = phenotype(lineID, traitID, traitval)
      insertedPhenoID = insert_phenotype(conn, pheno)
      insertedPhenoIDs.append(insertedPhenoID)
  return insertedPhenoIDs

def insert_phenotype(conn, phenotype):
  cur = conn.cursor()
  SQL = """INSERT INTO phenotype(phenotype_line, phenotype_trait, phenotype_value)
        VALUES (%s, %s, %s)
        ON CONFLICT DO NOTHING
        RETURNING phenotype_id;"""
  args_tuple = (phenotype.l, phenotype.t, phenotype.v)
  cur.execute(SQL, args_tuple)
  row = cur.fetchone()
  if row is not None:
    newID = row[0]
    conn.commit()
    cur.close()
    return newID
  else:
    return None

def insert_traits_from_traitlist(conn, traitlist):
  traitIDs = []
  for traitname in traitlist:
    traitObj = trait(traitname, None, None, None)
    insertedTraitID = insert_trait(conn, traitObj)
    traitIDs.append(insertedTraitID)
  return traitIDs

def insert_growout_type(conn, growout_type):
  cur = conn.cursor()
  SQL = """INSERT INTO growout_type(growout_type)
        VALUES (%s)
        ON CONFLICT DO NOTHING
        RETURNING growout_type_id;"""
  arg = (growout_type.t,)
  cur.execute(SQL, arg)
  row = cur.fetchone()
  if row is not None:
    newID = row[0]
    conn.commit()
    cur.close()
    return newID
  else:
    return None

def find_growout_type(conn, growout_type):
  cur = conn.cursor()
  cur.execute("SELECT growout_type_id FROM growout_type WHERE growout_type = '%s';" % growout_type)
  row = cur.fetchone()
  if row is not None:
    growout_type_ID = row[0]
    cur.close()
    return growout_type_ID
  else:
    return None

def insert_growout(conn, growout):
  cur = conn.cursor()
  SQL = """INSERT INTO growout(growout_name, growout_population, growout_location, year, growout_growout_type)
        VALUES (%s, %s, %s, %s, %s)
        ON CONFLICT DO NOTHING
        RETURNING growout_id;"""
  args_tuple = (growout.n, growout.p, growout.l, growout.y, growout.t)
  cur.execute(SQL, args_tuple)
  row = cur.fetchone()
  if row is not None:
    newID = row[0]
    conn.commit()
    cur.close()
    return newID
  else:
    return None

def insert_location(conn, location):
  cur = conn.cursor()
  SQL = """INSERT INTO location(country, state, city, code)
        VALUES (%s, %s, %s, %s)
        ON CONFLICT DO NOTHING
        RETURNING location_id;"""
  args_tuple = (location.c, location.s, location.i, location.o)
  cur.execute(SQL, args_tuple)
  row = cur.fetchone()
  if row is not None:
    newID = row[0]
    conn.commit()
    cur.close()
    return newID
  else:
    return None  
   
def find_location(conn, code):
  cur = conn.cursor()
  cur.execute("SELECT location_id FROM location WHERE code = '%s';" % code)
  row = cur.fetchone()
  if row is not None:
    location_ID = row[0]
    cur.close()
    return location_ID
  else:
    return None

def insert_gwas_algorithm(conn, gwas_algorithm):
  cur = conn.cursor()
  SQL = """INSERT INTO gwas_algorithm(gwas_algorithm)
        VALUES (%s)
        ON CONFLICT DO NOTHING
        RETURNING gwas_algorithm_id;"""
  args = (gwas_algorithm.a,)
  cur.execute(SQL, args)
  row = cur.fetchone()
  if row is not None:
    newID = row[0]
    conn.commit()
    cur.close()
    return newID
  else:
    return None

def insert_genotype_version(conn, genotype_version):
  cur = conn.cursor()
  SQL = """INSERT INTO genotype_version(genotype_version_name, genotype_version, reference_genome, genotype_version_population)
        VALUES (%s,%s,%s,%s)
        ON CONFLICT DO NOTHING
        RETURNING genotype_version_id;"""
  args_tuple = (genotype_version.n, genotype_version.v, genotype_version.r, genotype_version.p)
  cur.execute(SQL, args_tuple)
  row = cur.fetchone()
  if row is not None:
    newID = row[0]
    conn.commit()
    cur.close()
    return newID
  else:
    return None

def insert_imputation_method(conn, imputation_method):
  cur = conn.cursor()
  SQL = """INSERT INTO imputation_method(imputation_method)
        VALUES (%s)
        ON CONFLICT DO NOTHING
        RETURNING imputation_method_id;"""
  args = (imputation_method.m,)
  cur.execute(SQL, args)
  row = cur.fetchone()
  if row is not None:
    newID = row[0]
    conn.commit()
    cur.close()
    return newID
  else:
    return None

def insert_kinship_algorithm(conn, kinship_algorithm):
  cur = conn.cursor()
  SQL = """INSERT INTO kinship_algorithm(kinship_algorithm)
        VALUES (%s)
        ON CONFLICT DO NOTHING
        RETURNING kinship_algorithm_id;"""
  args = (kinship_algorithm.a,)
  cur.execute(SQL, args)
  row = cur.fetchone()
  if row is not None:
    newID = row[0]
    conn.commit()
    cur.close()
    return newID
  else:
    return None

def find_kinship_algorithm(conn, algorithm):
  cur = conn.cursor()
  cur.execute("SELECT kinship_algorithm_id FROM kinship_algorithm WHERE kinship_algorithm = '%s';" % algorithm)
  row = cur.fetchone()
  if row is not None:
    kinship_algorithm_ID = row[0]
    cur.close()
    return kinship_algorithm_ID
  else:
    return None

def insert_kinship(conn, kinship):
  cur = conn.cursor()
  SQL = """INSERT INTO kinship(kinship_algorithm, kinship_file_path)
        VALUES (%s,%s)
        ON CONFLICT DO NOTHING
        RETURNING kinship_id;"""
  args = (kinship.a, kinship.p)
  cur.execute(SQL, args)
  row = cur.fetchone()
  if row is not None:
    newID = row[0]
    conn.commit()
    cur.close()
    return newID
  else:
    return None

def insert_population_structure_algorithm(conn, population_structure_algorithm):
  cur = conn.cursor()
  SQL = """INSERT INTO population_structure_algorithm(population_structure_algorithm)
        VALUES (%s)
        ON CONFLICT DO NOTHING
        RETURNING population_structure_algorithm_id;"""
  args = (population_structure_algorithm.a,)
  cur.execute(SQL, args)
  row = cur.fetchone()
  if row is not None:
    newID = row[0]
    conn.commit()
    cur.close()
    return newID
  else:
    return None

def find_population_structure_algorithm(conn, algorithm):
  cur = conn.cursor()
  cur.execute("SELECT population_structure_algorithm_id FROM population_structure_algorithm WHERE population_structure_algorithm = '%s';" % algorithm)
  row = cur.fetchone()
  if row is not None:
    population_structure_algorithm_ID = row[0]
    cur.close()
    return population_structure_algorithm_ID
  else:
    return None

def insert_population_structure(conn, population_structure):
  cur = conn.cursor()
  SQL = """INSERT INTO population_structure(population_structure_algorithm, population_structure_file_path)
        VALUES (%s,%s)
        ON CONFLICT DO NOTHING
        RETURNING population_structure_id;"""
  args = (population_structure.a, population_structure.p)
  cur.execute(SQL, args)
  row = cur.fetchone()
  if row is not None:
    newID = row[0]
    conn.commit()
    cur.close()
    return newID
  else:
    return None

def find_gwas_algorithm(conn, gwas_algorithm):
  cur = conn.cursor()
  cur.execute("SELECT gwas_algorithm_id FROM gwas_algorithm WHERE gwas_algorithm = '%s';" % gwas_algorithm)
  row = cur.fetchone()
  if row is not None:
    gwas_algorithm_ID = row[0]
    cur.close()
    return gwas_algorithm_ID
  else:
    return None

def find_genotype_version(conn, genotype_version_name):
  cur = conn.cursor()
  cur.execute("SELECT genotype_version_id FROM genotype_version WHERE genotype_version_name = '%s';" % genotype_version_name)
  row = cur.fetchone()
  if row is not None:
    genotype_version_ID = row[0]
    cur.close()
    return genotype_version_ID
  else:
    return None

def find_imputation_method(conn, imputation_method):
  cur = conn.cursor()
  cur.execute("SELECT imputation_method_id FROM imputation_method WHERE imputation_method = '%s';" % imputation_method)
  row = cur.fetchone()
  if row is not None:
    imputation_method_ID = row[0]
    cur.close()
    return imputation_method_ID
  else:
    return None

def find_kinship(conn, kinship_file_path):
  cur = conn.cursor()
  cur.execute("SELECT kinship_id FROM kinship WHERE kinship_file_path = '%s';" % kinship_file_path)
  row = cur.fetchone()
  if row is not None:
    kinship_ID = row[0]
    cur.close()
    return kinship_ID
  else:
    return None

def find_population_structure(conn, population_structure_file_path):
  cur = conn.cursor()
  cur.execute("SELECT population_structure_id FROM population_structure WHERE population_structure_file_path = '%s';" % population_structure_file_path)
  row = cur.fetchone()
  if row is not None:
    population_structure_ID = row[0]
    cur.close()
    return population_structure_ID
  else:
    return None

def find_trait(conn, trait_name):
  cur = conn.cursor()
  cur.execute("SELECT trait_id FROM trait WHERE trait_name = '%s';" % trait_name)
  row = cur.fetchone()
  if row is not None:
    trait_ID = row[0]
    cur.close()
    return trait_ID
  else:
    return None

def parse_unique_runs_from_gwas_results_file(filepath):
  gwas_runs = []
  df = pd.read_csv(filepath)
  for index, row in df.iterrows():
    gwas_run = [row['trait'],row['nSNPs'],row['nLines']]
    if gwas_run not in gwas_runs:
      gwas_runs.append(gwas_run)
  return gwas_runs

def insert_gwas_run(conn, gwas_run):
  cur = conn.cursor()
  SQL = """INSERT INTO gwas_run(gwas_run_trait, nsnps, nlines, gwas_run_gwas_algorithm, gwas_run_genotype_version, missing_snp_cutoff_value, missing_line_cutoff_value, minor_allele_frequency_cutoff_value, gwas_run_imputation_method, gwas_run_kinship, gwas_run_population_structure)
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        ON CONFLICT DO NOTHING
        RETURNING gwas_run_id;"""
  args = (gwas_run.t, gwas_run.s, gwas_run.l, gwas_run.a, gwas_run.v, gwas_run.m, gwas_run.i, gwas_run.n, gwas_run.p, gwas_run.k, gwas_run.o)
  cur.execute(SQL, args)
  row = cur.fetchone()
  if row is not None:
    newID = row[0]
    conn.commit()
    cur.close()
    return newID
  else:
    return None

def insert_gwas_runs_from_gwas_results_file(conn, gwas_results_file, gwasRunAlgorithmID, gwasRunGenotypeVersionID, missing_snp_cutoff_value, missing_line_cutoff_value, minor_allele_frequency_cutoff_value, gwasRunImputationMethodID, gwasRunKinshipID, gwasRunPopulationStructureID):
  gwas_run_list = parse_unique_runs_from_gwas_results_file(gwas_results_file)
  insertedGwasRunIDs = []
  for gwas_run_item in gwas_run_list:
    traitID = find_trait(conn, gwas_run_item[0])
    gwas_run_obj = gwas_run(traitID, gwas_run_item[1], gwas_run_item[2], gwasRunAlgorithmID, gwasRunGenotypeVersionID, missing_snp_cutoff_value, missing_line_cutoff_value, minor_allele_frequency_cutoff_value, gwasRunImputationMethodID, gwasRunKinshipID, gwasRunPopulationStructureID)
    insertedGwasRunID = insert_gwas_run(conn, gwas_run_obj)
    insertedGwasRunIDs.append(insertedGwasRunID)
  return insertedGwasRunIDs

def find_gwas_run(conn, gwas_algorithm, missing_snp_cutoff_value, missing_line_cutoff_value, gwas_run_imputation_method, gwas_run_trait, nsnps, nlines, gwas_run_genotype_version, gwas_run_kinship, gwas_run_population_structure, minor_allele_frequency_cutoff_value):
  cur = conn.cursor()
  cur.execute("SELECT gwas_run_id FROM gwas_run WHERE gwas_run_gwas_algorithm = '%s' AND missing_snp_cutoff_value = '%s' AND missing_line_cutoff_value = '%s' AND gwas_run_imputation_method = '%s' AND gwas_run_trait = '%s' AND nsnps = '%s' AND nlines = '%s' AND gwas_run_genotype_version = '%s' AND gwas_run_kinship = '%s' AND gwas_run_population_structure = '%s' AND minor_allele_frequency_cutoff_value = '%s';" % (gwas_algorithm, missing_snp_cutoff_value, missing_line_cutoff_value, gwas_run_imputation_method, gwas_run_trait, nsnps, nlines, gwas_run_genotype_version, gwas_run_kinship, gwas_run_population_structure, minor_allele_frequency_cutoff_value))
  row = cur.fetchone()
  if row is not None:
    gwas_run_ID = row[0]
    cur.close()
    return gwas_run_ID
  else:
    return None

def insert_gwas_results_from_file(conn, speciesID, gwas_results_file, gwas_algorithm_ID, missing_snp_cutoff_value, missing_line_cutoff_value, imputationMethodID, genotypeVersionID, kinshipID, populationStructureID, minor_allele_frequency_cutoff_value):
  new_gwas_result_IDs = []
  df = pd.read_csv(gwas_results_file)
  for index, row in df.iterrows():
    trait = row['trait']
    traitID = find_trait(conn, trait)
    gwas_run_ID = find_gwas_run(conn, gwas_algorithm_ID, missing_snp_cutoff_value, missing_line_cutoff_value, imputationMethodID, traitID, row['nSNPs'], row['nLines'], genotypeVersionID, kinshipID, populationStructureID, minor_allele_frequency_cutoff_value)
    snp = row['SNP']
    snp_list = snp.split("_")
    chromosome = snp_list[0]
    chromosome = "chr"+str(chromosome)
    chromosomeID = find_chromosome(conn, chromosome, speciesID)
    basepair = snp_list[1]
    new_gwas_result = gwas_result(chromosomeID, basepair, gwas_run_ID, row['pval'], row['cofactor'], row['order'], row['nullPval'], row['modelAddedPval'], row['model'], row['PCs'])
    new_gwas_result_ID = insert_gwas_result(conn, new_gwas_result)
    new_gwas_result_IDs.append(new_gwas_result_ID)
  return new_gwas_result_IDs

def insert_gwas_result(conn, gwas_result):
  cur = conn.cursor()
  SQL = """INSERT INTO gwas_result(gwas_result_chromosome, basepair, gwas_result_gwas_run, pval, cofactor, _order, null_pval, model_added_pval, model, pcs)
        VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s, %s)
        ON CONFLICT DO NOTHING
        RETURNING gwas_result_id;"""
  args = (gwas_result.c, gwas_result.b, gwas_result.r, gwas_result.p, gwas_result.o, gwas_result.d, gwas_result.n, gwas_result.a, gwas_result.m, gwas_result.s)
  cur.execute(SQL, args)
  row = cur.fetchone()
  if row is not None:
    newID = row[0]
    conn.commit()
    cur.close()
    return newID
  else:
    return None

class species:
  def __init__(self, shortname, binomial, subspecies, variety):
    self.n = shortname
    self.b = binomial
    self.s = subspecies
    self.v = variety

class population:
  def __init__(self, population_name, population_species):
    self.n = population_name
    self.s = population_species

class line:
  def __init__(self, line_name, line_population):
    self.n = line_name
    self.p = line_population

class chromosome:
  def __init__(self, chromosome_name, chromosome_species):
    self.n = chromosome_name
    self.s = chromosome_species

class variant:
  def __init__(self, variant_species, variant_chromosome, variant_pos):
    self.s = variant_species
    self.c = variant_chromosome
    self.p = variant_pos

class genotype:
  def __init__(self, genotype_line, genotype_chromosome, genotype):
    self.l = genotype_line
    self.c = genotype_chromosome
    self.g = genotype

class trait:
  def __init__(self, trait_name, measurement_unit, measurement_device, description):
    self.n = trait_name
    self.u = measurement_unit
    self.m = measurement_device
    self.d = description

class phenotype:
  def __init__(self, phenotype_line, phenotype_trait, phenotype_value):
    self.l = phenotype_line
    self.t = phenotype_trait
    self.v = phenotype_value

class growout_type:
  def __init__(self, growout_type):
    self.t = growout_type

class growout:
  def __init__(self, growout_name, growout_population, growout_location, year, growout_growout_type):
    self.n = growout_name
    self.p = growout_population
    self.l = growout_location
    self.y = year
    self.t = growout_growout_type

class location:
  def __init__(self, country, state, city, code):
    self.c = country
    self.s = state
    self.i = city
    self.o = code

class gwas_algorithm:
  def __init__(self, gwas_algorithm):
    self.a = gwas_algorithm

class genotype_version:
  def __init__(self, genotype_version_name, genotype_version, reference_genome, genotype_version_population):
    self.n = genotype_version_name
    self.v = genotype_version
    self.r = reference_genome
    self.p = genotype_version_population

class imputation_method:
  def __init__(self, imputation_method):
    self.m = imputation_method

class kinship_algorithm:
  def __init__(self, kinship_algorithm):
    self.a = kinship_algorithm

class kinship:
  def __init__(self, kinship_algorithm, kinship_file_path):
    self.a = kinship_algorithm
    self.p = kinship_file_path

class population_structure_algorithm:
  def __init__(self, population_structure_algorithm):
    self.a = population_structure_algorithm

class population_structure:
  def __init__(self, population_structure_algorithm, population_structure_file_path):
    self.a = population_structure_algorithm
    self.p = population_structure_file_path

class gwas_run:
  def __init__(self, gwas_run_trait, nsnps, nlines, gwas_run_gwas_algorithm, gwas_run_genotype_version, missing_snp_cutoff_value, missing_line_cutoff_value, minor_allele_frequency_cutoff_value, gwas_run_imputation_method, gwas_run_kinship, gwas_run_population_structure):
    self.t = gwas_run_trait
    self.s = nsnps
    self.l = nlines
    self.a = gwas_run_gwas_algorithm
    self.v = gwas_run_genotype_version
    self.m = missing_snp_cutoff_value
    self.i = missing_line_cutoff_value
    self.n = minor_allele_frequency_cutoff_value
    self.p = gwas_run_imputation_method
    self.k = gwas_run_kinship
    self.o = gwas_run_population_structure

class gwas_result:
  def __init__(self, gwas_result_chromosome, basepair, gwas_result_gwas_run, pval, cofactor, _order, null_pval, model_added_pval, model, pcs):
    self.c = gwas_result_chromosome
    self.b = basepair
    self.r = gwas_result_gwas_run
    self.p = pval
    self.o = cofactor
    self.d = _order
    self.n = null_pval
    self.a = model_added_pval
    self.m = model
    self.s = pcs

if __name__ == '__main__':
  conn = connect()
  #########################################################
  # ADD A HARD-CODED SPECIES TO DB USING insert_species() #
  #########################################################
  #mySpecies = species('soybean', 'Glycine max', None, None)
  #insertedSpeciesID = insert_species(conn, mySpecies)
  #print(insertedSpeciesID)

  ###############################################################
  # ADD A HARD-CODED POPULATION TO DB USING insert_population() #
  ###############################################################
  #myPopulation = population('Maize282',maizeSpeciesID)
  #insertedPopulationID = insert_population(conn, myPopulation)
  #print(insertedPopulationID)

  ###########################################################
  # LOOK UP ID OF A HARD-CODED SPECIES USING find_species() #
  ###########################################################
  maizeSpeciesID = find_species(conn, 'maize')
  #print("SpeciesID of maize:")
  #print(maizeSpeciesID)

  #################################################################
  # LOOK UP ID OF A HARD-CODED POPULATION USING find_population() #
  #################################################################
  maize282popID = find_population(conn, 'Maize282')
  #print("PopulationID of Maize282:")
  #print(maize282popID)

  #################################################################
  # LOOK UP ID OF A HARD-CODED CHROMOSOME USING find_chromosome() #
  #################################################################
  #maizeChr10ID = find_chromosome(conn, 'chr10', maizeSpeciesID)
  #print("ChromosomeID of Maize Chr10:")
  #print(maizeChr10ID) 

  #####################################################
  # LOOK UP ID OF A HARD-CODED LINE USING find_line() #
  #####################################################
  B73lineID = find_line(conn, '282set_B73', maize282popID)
  
  ###################################################################
  # LOOK UP ID OF A HARD-CODED GROWOUT_TYPE USING find_chromosome() #
  ###################################################################
  fieldGrowoutTypeID = find_growout_type(conn, 'field')
  #print("fieldGrowoutTypeID:")
  #print(fieldGrowoutTypeID)
  
  ###############################################################
  # LOOK UP ID OF A HARD-CODED LOCATION USING find_chromosome() #
  ###############################################################
  PUlocID = find_location(conn, 'PU')
  NYlocID = find_location(conn, "NY")
  FLlocID = find_location(conn, "FL")
  PRlocID = find_location(conn, "PR")
  NClocID = find_location(conn, "NC")
  SAlocID = find_location(conn, "SA")
  MOlocID = find_location(conn, "MO")

  ########################################################
  # GET LINES FROM SPECIFIED 012.indv FILE AND ADD TO DB #
  ########################################################
  #insertedLineIDs = insert_lines_from_file(conn, '/home/mwohl/Downloads/GWASdata/chr10_282_agpv4.012.indv', maize282popID)
  #print("Inserted line IDs:")
  #print(insertedLineIDs)

  ###########################################
  # ADD ALL CHROMOSOMES FOR A SPECIES TO DB #
  ###########################################
  #insertedChromosomeIDs = insert_all_chromosomes_for_species(conn, 10, maizeSpeciesID)
  #print("Inserted chromosome IDs:")
  #print(insertedChromosomeIDs)

  ###########################################################
  # ADD ALL GENOTYPES FROM A ONE-CHROMOSOME .012 FILE TO DB #
  ###########################################################
  #insertedGenotypeIDs = insert_genotypes_from_file(conn,'/home/mwohl/Downloads/GWASdata/chr10_282_agpv4.012' , '/home/mwohl/Downloads/GWASdata/chr10_282_agpv4.012.indv', maizeChr10ID, maize282popID)
  #print("Inserted genotype IDs:")
  #print(insertedGenotypeIDs)


  ##################################################
  # GET VARIANTS FROM .012.pos FILE AND ADD TO  DB #
  ##################################################
  #insertedVariantIDs = insert_variants_from_file(conn, '/home/mwohl/Downloads/GWASdata/chr10_282_agpv4.012.pos', maizeSpeciesID, maizeChr10ID)
  #print("num inserted variants:")
  #print(len(insertedVariantIDs))

  ##################################################
  # PARSE TRAITS FROM PHENOTYPE FILE AND ADD TO DB #
  ##################################################
  #phenotypeRawData = pd.read_csv('/home/mwohl/Downloads/GWASdata/5.mergedWeightNorm.LM.rankAvg.longFormat.csv', index_col=0)
  #traits = list(phenotypeRawData)
  #insertedTraitIDs = insert_traits_from_traitlist(conn, traits)
  #print("num inserted traits:")
  #print(len(insertedTraitIDs))
  #print("Inserted trait IDs:")
  #print(insertedTraitIDs)
  
  ############################################
  # PARSE PHENOTYPES FROM FILE AND ADD TO DB #
  ############################################
  #insertedPhenoIDs = insert_phenotypes_from_file(conn, '/home/mwohl/Downloads/GWASdata/5.mergedWeightNorm.LM.rankAvg.longFormat.csv', maize282popID)
  #print("num phenotypes inserted:")
  #print(len(insertedPhenoIDs))
  #print("phenoIDs:")
  #print(insertedPhenoIDs)

  #########################################
  # ADD NEW HARD-CODED GROWOUT_TYPE TO DB #
  #########################################
  #greenhouse_GrowoutType = growout_type("greenhouse")
  #greenhouse_GrowoutTypeID = insert_growout_type(conn, greenhouse_GrowoutType)

  #phenotyper_GrowoutType = growout_type("phenotyper")
  #phenotyper_GrowoutTypeID = insert_growout_type(conn, phenotyper_GrowoutType)

  #field_GrowoutType = growout_type("field")
  #field_GrowoutTypeID = insert_growout_type(conn, field_GrowoutType)

  #####################################
  # ADD NEW HARD-CODED LOCATION TO DB #
  #####################################
  #newLocation = location("United States", "Indiana", "West Lafayette", "PU")
  #newLocation = location("United States", "New York", None, "NY")
  #newLocation = location("United States", "Florida", None, "FL")
  #newLocation = location("United States", "Puerto Rico", None, "PR")
  #newLocation = location("United States", "North Carolina", None, "NC")
  #newLocation = location("South Africa", None, None, "SA")
  #newLocation = location("United States", "Missouri", None, "MO")
  #newLocationID = insert_location(conn, newLocation)
  #print(newLocationID)

  ####################################
  # ADD NEW HARD-CODED GROWOUT TO DB #
  ####################################
  #newGrowout = growout("PU09", maize282popID, PUlocID, 2009, fieldGrowoutTypeID)
  #newGrowout = growout("NY06", maize282popID, NYlocID, 2006, fieldGrowoutTypeID)
  #newGrowout = growout("NY10", maize282popID, NYlocID, 2010, fieldGrowoutTypeID)
  #newGrowout = growout("FL06", maize282popID, FLlocID, 2006, fieldGrowoutTypeID)
  #newGrowout = growout("PR06", maize282popID, PRlocID, 2006, fieldGrowoutTypeID)
  #newGrowout = growout("NC06", maize282popID, NClocID, 2006, fieldGrowoutTypeID)
  #newGrowout = growout("PU10", maize282popID, PUlocID, 2010, fieldGrowoutTypeID)
  #newGrowout = growout("SA06", maize282popID, SAlocID, 2006, fieldGrowoutTypeID)
  #newGrowout = growout("MO06", maize282popID, MOlocID, 2006, fieldGrowoutTypeID)
  #newGrowoutID = insert_growout(conn, newGrowout)
  #print("new growout's ID:")
  #print(newGrowoutID)
  
  ###########################################
  # ADD NEW HARD-CODED GWAS_ALGORITHM TO DB #
  ###########################################
  #newGWASalgorithm = gwas_algorithm("MLMM")
  #newGWASalgorithm = gwas_algorithm("EMMAx")
  #newGWASalgorithm = gwas_algorithm("GAPIT")
  #newGWASalgorithm = gwas_algorithm("FarmCPU")
  #newGWASalgorithmID = insert_gwas_algorithm(conn, newGWASalgorithm)
  #print("new GWAS algorithm ID:")
  #print(newGWASalgorithmID)

  #############################################
  # ADD NEW HARD-CODED GENOTYPE_VERSION TO DB #
  #############################################
  #newGenotypeVersion = genotype_version("maize282_agpv4_B73", "apgv4", B73lineID, maize282popID)
  #newGenotypeVersionID = insert_genotype_version(conn, newGenotypeVersion)
  #print("New genotype_version ID:")
  #print(newGenotypeVersionID)

  ##############################################
  # ADD NEW HARD-CODED IMPUTATION_METHOD TO DB #
  ##############################################
  #newImputationMethod = imputation_method("impute to major allele")
  #newImputationMethod = imputation_method("impute to minor allele")
  #newImputationMethod = imputation_method("impute to average allele")
  #newImputationMethod = imputation_method("IMPUTE")
  #newImputationMethod = imputation_method("BEAGLE")
  #newImputationMethodID = insert_imputation_method(conn, newImputationMethod)
  #print("Imputatin Method ID:")
  #print(newImputationMethodID)

  ##############################################
  # ADD NEW HARD-CODED KINSHIP_ALGORITHM TO DB #
  ##############################################
  #newKinshipAlgorithm = kinship_algorithm("loiselle")
  #newKinshipAlgorithm = kinship_algorithm("van raden")
  #newKinshipAlgorithm = kinship_algorithm("Synbreed_realizedAB")
  #newKinshipAlgorithmID = insert_kinship_algorithm(conn, newKinshipAlgorithm)
  #print("Kinship Algorithm ID:")
  #print(newKinshipAlgorithmID)

  ###############################################################################
  # LOOK UP ID OF A HARD-CODED KINSHIP_ALGORITHM USING find_kinship_algorithm() #
  ###############################################################################
  #VanRadenID = find_kinship_algorithm(conn, "van raden")
  #print("Van Raden kinship alg ID:")
  #print(VanRadenID)  

  ####################################
  # ADD NEW HARD-CODED KINSHIP TO DB #
  ####################################
  #newKinship = kinship(VanRadenID, "/opt/BaxDB/file_storage/kinship_files/4.AstleBalding.synbreed.kinship.csv")
  #newKinshipID = insert_kinship(conn, newKinship)
  #print("New kinship ID:")
  #print(newKinshipID)

  ###########################################################
  # ADD NEW HARD-CODED POPULATION_STRUCTURE_ALGORITHM TO DB #
  ###########################################################
  #newPopulationStructureAlgorithm = population_structure_algorithm("Eigenstrat")
  #newPopulationStructureAlgorithm = population_structure_algorithm("STRUCTURE")
  #newPopulationStructureAlgorithm = population_structure_algorithm("FastSTRUCTURE")
  #newPopulationStructureAlgorithmID = insert_population_structure_algorithm(conn, newPopulationStructureAlgorithm)
  #print("pop structure algorithm ID:")
  #print(newPopulationStructureAlgorithmID)

  #########################################################################################################
  # LOOK UP ID OF A HARD-CODED POPULATION_STRUCTURE_ALGORITHM USING find_population_structure_algorithm() #
  #########################################################################################################
  #EigenstratID = find_population_structure_algorithm(conn, "Eigenstrat")
  #print("Eigenstrat pop str alg ID:")
  #print(EigenstratID)

  #################################################
  # ADD NEW HARD-CODED POPULATION_STRUCTURE TO DB #
  #################################################
  #newPopulationStructure = population_structure(EigenstratID, "/opt/BaxDB/file_storage/population_structure_files/4.Eigenstrat.population.structure.10PCs.csv")
  #newPopulationStructureID = insert_population_structure(conn, newPopulationStructure)
  #print("New population structure ID:")
  #print(newPopulationStructureID)

  #############################################
  # LOOK UP ID OF A HARD-CODED GWAS_ALGORITHM #
  #############################################
  MLMMalgorithmID = find_gwas_algorithm(conn, "MLMM")
  #print("MLMM algorithm ID:")
  #print(MLMMalgorithmID)

  ###############################################
  # LOOK UP ID OF A HARD-CODED GENOTYPE_VERSION #
  ###############################################
  B73_agpv4_maize282_versionID = find_genotype_version(conn, "B73 RefGen_v4_AGPv4_Maize282")
  #print("B73 agpv4 maize282 genotype version: ")
  #print(B73_agpv4_maize282_versionID)  

  ################################################
  # LOOK UP ID OF A HARD-CODED IMPUTATION_METHOD #
  ################################################
  majorAlleleImputationID = find_imputation_method(conn, "impute to major allele")
  #print("major allele imputation ID: ")
  #print(majorAlleleImputationID)  

  ######################################
  # LOOK UP ID OF A HARD-CODED KINSHIP #
  ######################################
  kinshipID = find_kinship(conn, "/opt/BaxDB/file_storage/kinship_files/4.AstleBalding.synbreed.kinship.csv")
  #print("kinshipID: ")
  #print(kinshipID)  

  ###################################################
  # LOOK UP ID OF A HARD-CODED POPULATION_STRUCTURE #
  ###################################################
  populationStructureID = find_population_structure(conn, "/opt/BaxDB/file_storage/population_structure_files/4.Eigenstrat.population.structure.10PCs.csv")
  #print("population structure ID: ")
  #print(populationStructureID)

  ###########################################
  # PARSE GWAS_RUNS FROM FILE AND ADD TO DB #
  ###########################################
  #insertedGwasRunIDs = insert_gwas_runs_from_gwas_results_file(conn, '/home/mwohl/Downloads/GWASdata/9.mlmmResults.csv', MLMMalgorithmID, B73_agpv4_maize282_versionID, 0.2, 0.2, 0.1, majorAlleleImputationID, kinshipID, populationStructureID)
  #print("Inserted gwas_run IDs:")
  #print(insertedGwasRunIDs)

  #######################################
  # LOOK UP ID OF A HARD-CODED GWAS_RUN #
  #######################################
   #gwasRunID = find_gwas_run(conn, MLMMalgorithmID, 0.2, 0.2, majorAlleleImputationID, ---, ---, ---, B73_agpv4_maize282_versionID, kinshipID, populationStructureID, 0.1) 

  ##############################################
  # PARSE GWAS_RESULTS FROM FILE AND ADD TO DB #
  ##############################################
  insertedGwasResultIDs = insert_gwas_results_from_file(conn, maizeSpeciesID, '/home/mwohl/Downloads/GWASdata/9.mlmmResults.csv', MLMMalgorithmID, 0.2, 0.2, majorAlleleImputationID, B73_agpv4_maize282_versionID, kinshipID, populationStructureID, 0.1)
  #print("Inserted gwas result IDs: ")
  #print(insertedGwasResultIDs)
