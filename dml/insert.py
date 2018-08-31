# Insert objects (as defined in models.py) into the database

import pandas as pd
import numpy as np
import parsinghelpers as ph
from models import species, population, line, chromosome, variant, genotype, trait, phenotype, growout_type, growout, location, gwas_algorithm, genotype_version, imputation_method, kinship_algorithm, kinship, population_structure_algorithm, population_structure, gwas_run, gwas_result

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
  chrlist = ph.generate_chromosome_list(numChromosomes)
  insertedChromosomeIDs = []
  for chrname in chrlist:
    chrobj = chromosome(chrname, speciesID)
    insertedChromosomeID = insert_chromosome(conn, chrobj)
    insertedChromosomeIDs.append(insertedChromosomeID)
  return insertedChromosomeIDs

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
  linelist = ph.parse_lines_from_file(lineFile)
  insertedLineIDs = []
  for linename in linelist:
    lineobj = line(linename, populationID)
    insertedLineID = insert_line(conn, lineobj)
    insertedLineIDs.append(insertedLineID)
  return insertedLineIDs

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
  variantlist = ph.parse_variants_from_file(variantPosFile)
  print('num variants:')
  print(len(variantlist))
  insertedVariantIDs = []
  for variantpos in variantlist:
    variantobj = variant(speciesID, chromosomeID, variantpos)
    insertedVariantID = insert_variant(conn, variantobj)
    insertedVariantIDs.append(insertedVariantID)
  return insertedVariantIDs

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

def insert_genotypes_from_file(conn, genotypeFile, lineFile, chromosomeID, populationID):
  genotypes = ph.parse_genotypes_from_file(genotypeFile)
  linelist = ph.parse_lines_from_file(lineFile)
  lineIDlist = ph.convert_linelist_to_lineIDlist(conn, linelist, populationID)
  zipped = zip(lineIDlist, genotypes)
  ziplist = list(zipped)
  insertedGenotypeIDs = []
  for zippedpair in ziplist:
    genotypeObj = genotype(zippedpair[0], chromosomeID, zippedpair[1])
    insertedGenotypeID = insert_genotype(conn, genotypeObj)
    insertedGenotypeIDs.append(insertedGenotypeID)
  return insertedGenotypeIDs

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
  gwas_run_list = ph.parse_unique_runs_from_gwas_results_file(gwas_results_file)
  insertedGwasRunIDs = []
  for gwas_run_item in gwas_run_list:
    traitID = find_trait(conn, gwas_run_item[0])
    gwas_run_obj = gwas_run(traitID, gwas_run_item[1], gwas_run_item[2], gwasRunAlgorithmID, gwasRunGenotypeVersionID, missing_snp_cutoff_value, missing_line_cutoff_value, minor_allele_frequency_cutoff_value, gwasRunImputationMethodID, gwasRunKinshipID, gwasRunPopulationStructureID)
    insertedGwasRunID = insert_gwas_run(conn, gwas_run_obj)
    insertedGwasRunIDs.append(insertedGwasRunID)
  return insertedGwasRunIDs

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
    
    pcs = row['PCs']
    if type(pcs) == str:
      pcs_list = pcs.split(":")
      pcs_list = [int(x) for x in pcs_list]
    elif np.isnan(pcs):
      pcs_list = None
 
    new_gwas_result = gwas_result(chromosomeID, basepair, gwas_run_ID, row['pval'], row['cofactor'], row['order'], row['nullPval'], row['modelAddedPval'], row['model'], pcs_list)
    if new_gwas_result:
      print("yep")
    else:
      print("nope")
    new_gwas_result_ID = insert_gwas_result(conn, new_gwas_result)
    new_gwas_result_IDs.append(new_gwas_result_ID)
  return new_gwas_result_IDs