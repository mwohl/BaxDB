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

if __name__ == '__main__':
  conn = connect()
  #mySpecies = species('maize', 'Zea mays', None, None)
  #insertedSpeciesID = insert_species(conn, mySpecies)
  #print(insertedSpeciesID)

  maizeSpeciesID = find_species(conn, 'maize')
  print(maizeSpeciesID)

  #myPopulation = population('Maize282',maizeSpeciesID)
  #insertedPopulationID = insert_population(conn, myPopulation)
  #print(insertedPopulationID)

  chrlist = []
  for count in range(1,11):
    chrname = 'chr'+str(count)
    chrom = chromosome(chrname,maizeSpeciesID)
    chrlist.append(chrom)
  print(chrlist)
  for item in chrlist:
    insertedChromosomeID = insert_chromosome(conn, item)
    print(insertedChromosomeID)
