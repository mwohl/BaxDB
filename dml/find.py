# Look up IDs of various elements in the database

import pandas as pd
import numpy as np
import psycopg2
import csv
import insert
from dbconnect import config, connect
from models import species, population, line, chromosome, variant, genotype, trait, phenotype, growout_type, growout, location, gwas_algorithm, genotype_version, imputation_method, kinship_algorithm, kinship, population_structure_algorithm, population_structure, gwas_run, gwas_result

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