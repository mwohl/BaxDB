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

def generate_chromosome_list(numChromosomes):
  chrlist = []
  for count in range(1,numChromosomes+1):
    chrname = 'chr'+str(count)
    chrlist.append(chrname)
  return chrlist

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

def parse_unique_runs_from_gwas_results_file(filepath):
  gwas_runs = []
  df = pd.read_csv(filepath)
  for index, row in df.iterrows():
    gwas_run = [row['trait'],row['nSNPs'],row['nLines']]
    if gwas_run not in gwas_runs:
      gwas_runs.append(gwas_run)
  return gwas_runs

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

if __name__ == '__main__':
  conn = connect()
  #########################################################
  # ADD A HARD-CODED SPECIES TO DB USING insert_species() #
  #########################################################
  mySpecies = species('testSpecies', 'Test sp', None, None)
  insertedSpeciesID = insert.insert_species(conn, mySpecies)
  print(insertedSpeciesID)

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
  #insertedGwasResultIDs = insert_gwas_results_from_file(conn, maizeSpeciesID, '/home/mwohl/Downloads/GWASdata/9.mlmmResults.csv', MLMMalgorithmID, 0.2, 0.2, majorAlleleImputationID, B73_agpv4_maize282_versionID, kinshipID, populationStructureID, 0.1)
  #print("Inserted gwas result IDs: ")
  #print(insertedGwasResultIDs)
