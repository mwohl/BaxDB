from configparser import ConfigParser
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
        VALUES (%s, %s, %s, %s) RETURNING species_id;"""
  args_tuple = (species.n, species.b, species.s, species.v)
  cur.execute(SQL, args_tuple)
  newID = cur.fetchone()[0]
  conn.commit()
  cur.close()
  return newID

def find_species(conn, speciesShortname):
  cur = conn.cursor()
  cur.execute("SELECT species_id FROM species WHERE shortname = '%s';" %  speciesShortname)
  speciesID = cur.fetchone()[0]
  cur.close()
  return speciesID

def insert_population(conn, population): 
  cur = conn.cursor()
  SQL = """INSERT INTO population (population_name, population_species)
        VALUES (%s, %s) RETURNING population_id;"""
  args_tuple = (population.n, population.s)
  cur.execute(SQL, args_tuple)
  newID = cur.fetchone()[0]
  conn.commit()
  cur.close()
  return newID  

def find_population(conn, populationName):
  cur = conn.cursor()
  cur.execute("SELECT population_id FROM population WHERE population_name = '%s';" % populationName)
  populationID = cur.fetchone()[0]
  cur.close()
  return populationID

def generate_chromosome_list(numChromosomes):
  chrlist = []
  for count in range(1,numChromosomes+1):
    chrname = 'chr'+str(count)
    chrlist.append(chrname)
  return chrlist

def insert_chromosome(conn, chromosome):
  cur = conn.cursor()
  SQL = """INSERT INTO chromosome (chromosome_name, chromosome_species)
        VALUES (%s, %s) RETURNING chromosome_id;"""
  args_tuple = (chromosome.n, chromosome.s)
  cur.execute(SQL, args_tuple)
  newID = cur.fetchone()[0]
  conn.commit()
  cur.close()
  return newID

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
  chromosomeID = cur.fetchone()[0]
  cur.close()
  return chromosomeID

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
        VALUES (%s, %s) RETURNING line_id;"""
  args_tuple = (line.n, line.p)
  cur.execute(SQL, args_tuple)
  newID = cur.fetchone()[0]
  conn.commit()
  cur.close()
  return newID

def insert_lines_from_file(conn, lineFile, populationID):
  linelist = parse_lines_from_file(lineFile)
  insertedLineIDs = []
  for linename in linelist:
    lineobj = line(linename, populationID)
    insertedLineID = insert_line(conn, lineobj)
    insertedLineIDs.append(insertedLineID)
  return insertedLineIDs

def find_line(conn, line_name):
  cur = conn.cursor()
  cur.execute("SELECT line_id FROM line WHERE line_name = '%s';" % line_name)
  lineID = cur.fetchone()[0]
  cur.close()
  return lineID

def convert_linelist_to_lineIDlist(conn, linelist):
  lineIDlist = []
  for linename in linelist:
    lineID = find_line(conn, linename)
    lineIDlist.append(lineID)
  return lineIDlist

def parse_variants_from_file(variantPosFile):
  with open(variantPosFile) as f:
    variantReader = csv.reader(f, delimiter='\t')
    variantlist = []
    for variant in variantReader:
      variantlist.append(variant[1])
  return variantlist

#NEEDS TESTING
def insert_variant(conn, variant):
  cur = conn.cursor()
  SQL = """INSERT INTO variant(variant_species, variant_chromosome, variant_pos)
        VALUES (%s,%s,%s) RETURNING variant_id;"""
  args_tuple = (variant.s, variant.c, variant.p)
  cur.execute(SQL, args_tuple)
  newID = cur.fetchone()[0]
  conn.commit()
  cur.close()
  return newID

#NEEDS TESTING
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

def insert_genotypes_from_file(conn, genotypeFile, lineFile, chromosomeID):
  genotypes = parse_genotypes_from_file(genotypeFile)
  linelist = parse_lines_from_file(lineFile)
  lineIDlist = convert_linelist_to_lineIDlist(conn, linelist)
  zipped = zip(lineIDlist, genotypes)
  ziplist = list(zipped)
  insertedGenotypeIDs = []
  for zippedpair in ziplist:
    print(type(zippedpair))
    print(len(zippedpair))
    print(zippedpair)
    genotypeObj = genotype(zippedpair[0], chromosomeID, zippedpair[1])
    insertedGenotypeID = insert_genotype(conn, genotypeObj)
    insertedGenotypeIDs.append(insertedGenotypeID)
  return insertedGenotypeIDs

def insert_genotype(conn, genotype):
  cur = conn.cursor()
  SQL = """INSERT INTO genotype(line_ref, chromosome_id, genotype)
        VALUES (%s,%s,%s) RETURNING line_ref;"""
  args_tuple = (genotype.l, genotype.c, genotype.g)
  cur.execute(SQL, args_tuple)
  newID = cur.fetchone()[0]
  conn.commit()
  cur.close()
  return newID

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
  def __init__(self, line_ref, chromosome_id, genotype):
    self.l = line_ref
    self.c = chromosome_id
    self.g = genotype

if __name__ == '__main__':
  conn = connect()
  #########################################################
  # ADD A HARD-CODED SPECIES TO DB USING insert_species() #
  #########################################################
  #mySpecies = species('maize', 'Zea mays', None, None)
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
  print("SpeciesID of maize:")
  print(maizeSpeciesID)

  #################################################################
  # LOOK UP ID OF A HARD-CODED POPULATION USING find_population() #
  #################################################################
  maize282popID = find_population(conn, 'Maize282')
  print("PopulationID of Maize282:")
  print(maize282popID)

  #################################################################
  # LOOK UP ID OF A HARD-CODED CHROMOSOME USING find_chromosome() #
  #################################################################
  maizeChr2ID = find_chromosome(conn, 'chr2', maizeSpeciesID)
  print("ChromosomeID of Maize Chr2:")
  print(maizeChr2ID) 

  ########################################################
  # GET LINES FROM SPECIFIED 012.indv FILE AND ADD TO DB #
  ########################################################
  #insertedLineIDs = insert_lines_from_file(conn, '/home/mwohl/Downloads/GWASdata/chr2_282_agpv4.012.indv', maize282popID)
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
  insertedGenotypeIDs = insert_genotypes_from_file(conn,'/home/mwohl/Downloads/GWASdata/chr2_282_agpv4.012' , '/home/mwohl/Downloads/GWASdata/chr2_282_agpv4.012.indv', maizeChr2ID)
  print("Inserted genotype IDs:")
  print(insertedGenotypeIDs)


  ######################################################################################
  # IN PROGRESS... code for getting variants from .012.pos files and inserting into DB #
  ######################################################################################

  # NEEDS TESTING (make sure to change above chromosome calling to chr1)
  #insertedVariantIDs = insert_variants_from_file(conn, '/home/mwohl/Downloads/GWASdata/chr1_282_agpv4.012.pos', maizeSpeciesID, maizeChr1ID)
  #print("num inserted variants:")
  #print(len(insertedVariantIDs))
