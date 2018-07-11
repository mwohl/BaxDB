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
  arg = (growout_type.t)
  cur.execute(SQL, arg)
  row = cur.fetchone()
  if row is not None:
    newID = row[0]
    conn.commit()
    cur.close()
    return newID
  else:
    return None

def insert_growout(conn, growout):
  cur = conn.cursor()
  SQL = """INSERT INTO growout(growout_name, growout_population, growout_location, growout_year, growout_type)
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
  def __init__(self, growout_name, growout_population, growout_location, growout_year, growout_type)
    self.n = growout_name
    self.p = growout_population
    self.l = growout_location
    self.y = growout_year
    self.t = growout_type

class location:
  def __init__(self, location_country, location_state, loctaion_city):
    self.c = location_country
    self.s = location_state
    self.i = location_city

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
  #maizeSpeciesID = find_species(conn, 'maize')
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
  #Mo17_lineID = find_line(conn, '282set_Mo17', maize282popID)
  #print("LineID of Mo17:")
  #print(Mo17_lineID)

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

  ####################################
  # ADD NEW HARD-CODED GROWOUT TO DB #
  ####################################
  #newGrowout = growout()
