# This script contains most of the code executed to insert all of the Maize282 dataset into the GWAS database.
# It uses functions from the other modules in the BaxDB code -- insert functions from insert.py, find functions from find.py, helper functions from parsinghelpers.py, classes defined in models.py, and connection/configuration functions from dbconnect.py

import pandas as pd
import numpy as np
import psycopg2
import csv
import insert
import find
from dbconnect import config, connect
from models import species, population, line, chromosome, variant, genotype, trait, phenotype, growout_type, growout, location, gwas_algorithm, genotype_version, imputation_method, kinship_algorithm, kinship, population_structure_algorithm, population_structure, gwas_run, gwas_result


if __name__ == '__main__':
  conn = connect()

  # ADD A HARD-CODED SPECIES TO DB USING insert_species()
  mySpecies = species('maize', 'Zea mays', None, None)
  insertedSpeciesID = insert.insert_species(conn, mySpecies)
  print(insertedSpeciesID)

  # ADD A HARD-CODED POPULATION TO DB USING insert_population()
  myPopulation = population('Maize282',maizeSpeciesID)
  insertedPopulationID = insert.insert_population(conn, myPopulation)
  print(insertedPopulationID)

  # LOOK UP ID OF A HARD-CODED SPECIES USING find_species()
  maizeSpeciesID = find.find_species(conn, 'maize')
  print("SpeciesID of maize:")
  print(maizeSpeciesID)

  # LOOK UP ID OF A HARD-CODED POPULATION USING find_population()
  maize282popID = find.find_population(conn, 'Maize282')
  print("PopulationID of Maize282:")
  print(maize282popID)

  # LOOK UP ID OF A HARD-CODED CHROMOSOME USING find_chromosome()
  maizeChr10ID = find.find_chromosome(conn, 'chr10', maizeSpeciesID)
  print("ChromosomeID of Maize Chr10:")
  print(maizeChr10ID) 

  # LOOK UP ID OF A HARD-CODED LINE USING find_line()
  B73lineID = find.find_line(conn, '282set_B73', maize282popID)
  
  # LOOK UP ID OF A HARD-CODED GROWOUT_TYPE USING find_chromosome()
  fieldGrowoutTypeID = find.find_growout_type(conn, 'field')
  print("fieldGrowoutTypeID:")
  print(fieldGrowoutTypeID)
  
  # LOOK UP ID OF A HARD-CODED LOCATION USING find_chromosome()
  PUlocID = find.find_location(conn, 'PU')
  NYlocID = find.find_location(conn, "NY")
  FLlocID = find.find_location(conn, "FL")
  PRlocID = find.find_location(conn, "PR")
  NClocID = find.find_location(conn, "NC")
  SAlocID = find.find_location(conn, "SA")
  MOlocID = find.find_location(conn, "MO")

  # GET LINES FROM SPECIFIED 012.indv FILE AND ADD TO DB
  insertedLineIDs = insert.insert_lines_from_file(conn, '/home/mwohl/Downloads/GWASdata/chr10_282_agpv4.012.indv', maize282popID)
  print("Inserted line IDs:")
  print(insertedLineIDs)

  # ADD ALL CHROMOSOMES FOR A SPECIES TO DB
  insertedChromosomeIDs = insert.insert_all_chromosomes_for_species(conn, 10, maizeSpeciesID)
  print("Inserted chromosome IDs:")
  print(insertedChromosomeIDs)

  # ADD ALL GENOTYPES FROM A ONE-CHROMOSOME .012 FILE TO DB
  insertedGenotypeIDs = insert.insert_genotypes_from_file(conn,'/home/mwohl/Downloads/GWASdata/chr10_282_agpv4.012' , '/home/mwohl/Downloads/GWASdata/chr10_282_agpv4.012.indv', maizeChr10ID, maize282popID)
  print("Inserted genotype IDs:")
  print(insertedGenotypeIDs)

  # GET VARIANTS FROM .012.pos FILE AND ADD TO  DB
  insertedVariantIDs = insert.insert_variants_from_file(conn, '/home/mwohl/Downloads/GWASdata/chr10_282_agpv4.012.pos', maizeSpeciesID, maizeChr10ID)
  print("num inserted variants:")
  print(len(insertedVariantIDs))

  # PARSE TRAITS FROM PHENOTYPE FILE AND ADD TO DB
  phenotypeRawData = pd.read_csv('/home/mwohl/Downloads/GWASdata/5.mergedWeightNorm.LM.rankAvg.longFormat.csv', index_col=0)
  traits = list(phenotypeRawData)
  insertedTraitIDs = insert.insert_traits_from_traitlist(conn, traits)
  print("num inserted traits:")
  print(len(insertedTraitIDs))
  print("Inserted trait IDs:")
  print(insertedTraitIDs)
  
  # PARSE PHENOTYPES FROM FILE AND ADD TO DB
  insertedPhenoIDs = insert.insert_phenotypes_from_file(conn, '/home/mwohl/Downloads/GWASdata/5.mergedWeightNorm.LM.rankAvg.longFormat.csv', maize282popID)
  print("num phenotypes inserted:")
  print(len(insertedPhenoIDs))
  print("phenoIDs:")
  print(insertedPhenoIDs)

  # ADD NEW HARD-CODED GROWOUT_TYPE TO DB
  greenhouse_GrowoutType = growout_type("greenhouse")
  greenhouse_GrowoutTypeID = insert.insert_growout_type(conn, greenhouse_GrowoutType)

  phenotyper_GrowoutType = growout_type("phenotyper")
  phenotyper_GrowoutTypeID = insert.insert_growout_type(conn, phenotyper_GrowoutType)

  field_GrowoutType = growout_type("field")
  field_GrowoutTypeID = insert.insert_growout_type(conn, field_GrowoutType)

  # ADD NEW HARD-CODED LOCATION TO DB
  newLocation = location("United States", "Indiana", "West Lafayette", "PU")
  #newLocation = location("United States", "New York", None, "NY")
  #newLocation = location("United States", "Florida", None, "FL")
  #newLocation = location("United States", "Puerto Rico", None, "PR")
  #newLocation = location("United States", "North Carolina", None, "NC")
  #newLocation = location("South Africa", None, None, "SA")
  #newLocation = location("United States", "Missouri", None, "MO")
  newLocationID = insert.insert_location(conn, newLocation)
  print(newLocationID)

  # ADD NEW HARD-CODED GROWOUT TO DB
  newGrowout = growout("PU09", maize282popID, PUlocID, 2009, fieldGrowoutTypeID)
  #newGrowout = growout("NY06", maize282popID, NYlocID, 2006, fieldGrowoutTypeID)
  #newGrowout = growout("NY10", maize282popID, NYlocID, 2010, fieldGrowoutTypeID)
  #newGrowout = growout("FL06", maize282popID, FLlocID, 2006, fieldGrowoutTypeID)
  #newGrowout = growout("PR06", maize282popID, PRlocID, 2006, fieldGrowoutTypeID)
  #newGrowout = growout("NC06", maize282popID, NClocID, 2006, fieldGrowoutTypeID)
  #newGrowout = growout("PU10", maize282popID, PUlocID, 2010, fieldGrowoutTypeID)
  #newGrowout = growout("SA06", maize282popID, SAlocID, 2006, fieldGrowoutTypeID)
  #newGrowout = growout("MO06", maize282popID, MOlocID, 2006, fieldGrowoutTypeID)
  newGrowoutID = insert.insert_growout(conn, newGrowout)
  print("new growout's ID:")
  print(newGrowoutID)
  
  # ADD NEW HARD-CODED GWAS_ALGORITHM TO DB
  newGWASalgorithm = gwas_algorithm("MLMM")
  #newGWASalgorithm = gwas_algorithm("EMMAx")
  #newGWASalgorithm = gwas_algorithm("GAPIT")
  #newGWASalgorithm = gwas_algorithm("FarmCPU")
  newGWASalgorithmID = insert.insert_gwas_algorithm(conn, newGWASalgorithm)
  print("new GWAS algorithm ID:")
  print(newGWASalgorithmID)

  # ADD NEW HARD-CODED GENOTYPE_VERSION TO DB
  newGenotypeVersion = genotype_version("maize282_agpv4_B73", "apgv4", B73lineID, maize282popID)
  newGenotypeVersionID = insert.insert_genotype_version(conn, newGenotypeVersion)
  print("New genotype_version ID:")
  print(newGenotypeVersionID)

  # ADD NEW HARD-CODED IMPUTATION_METHOD TO DB
  newImputationMethod = imputation_method("impute to major allele")
  #newImputationMethod = imputation_method("impute to minor allele")
  #newImputationMethod = imputation_method("impute to average allele")
  #newImputationMethod = imputation_method("IMPUTE")
  #newImputationMethod = imputation_method("BEAGLE")
  newImputationMethodID = insert.insert_imputation_method(conn, newImputationMethod)
  print("Imputatin Method ID:")
  print(newImputationMethodID)

  
  # ADD NEW HARD-CODED KINSHIP_ALGORITHM TO DB
  newKinshipAlgorithm = kinship_algorithm("loiselle")
  #newKinshipAlgorithm = kinship_algorithm("van raden")
  #newKinshipAlgorithm = kinship_algorithm("Synbreed_realizedAB")
  newKinshipAlgorithmID = insert.insert_kinship_algorithm(conn, newKinshipAlgorithm)
  print("Kinship Algorithm ID:")
  print(newKinshipAlgorithmID)

  # LOOK UP ID OF A HARD-CODED KINSHIP_ALGORITHM USING find_kinship_algorithm()
  VanRadenID = find.find_kinship_algorithm(conn, "van raden")
  print("Van Raden kinship alg ID:")
  print(VanRadenID)  

  # ADD NEW HARD-CODED KINSHIP TO DB
  newKinship = kinship(VanRadenID, "/opt/BaxDB/file_storage/kinship_files/4.AstleBalding.synbreed.kinship.csv")
  newKinshipID = insert.insert_kinship(conn, newKinship)
  print("New kinship ID:")
  print(newKinshipID)

  # ADD NEW HARD-CODED POPULATION_STRUCTURE_ALGORITHM TO DB
  newPopulationStructureAlgorithm = population_structure_algorithm("Eigenstrat")
  #newPopulationStructureAlgorithm = population_structure_algorithm("STRUCTURE")
  #newPopulationStructureAlgorithm = population_structure_algorithm("FastSTRUCTURE")
  newPopulationStructureAlgorithmID = insert.insert_population_structure_algorithm(conn, newPopulationStructureAlgorithm)
  print("pop structure algorithm ID:")
  print(newPopulationStructureAlgorithmID)

  # LOOK UP ID OF A HARD-CODED POPULATION_STRUCTURE_ALGORITHM USING find_population_structure_algorithm()
  EigenstratID = find.find_population_structure_algorithm(conn, "Eigenstrat")
  print("Eigenstrat pop str alg ID:")
  print(EigenstratID)

  # ADD NEW HARD-CODED POPULATION_STRUCTURE TO DB
  newPopulationStructure = population_structure(EigenstratID, "/opt/BaxDB/file_storage/population_structure_files/4.Eigenstrat.population.structure.10PCs.csv")
  newPopulationStructureID = insert.insert_population_structure(conn, newPopulationStructure)
  print("New population structure ID:")
  print(newPopulationStructureID)

  # LOOK UP ID OF A HARD-CODED GWAS_ALGORITHM
  MLMMalgorithmID = find.find_gwas_algorithm(conn, "MLMM")
  print("MLMM algorithm ID:")
  print(MLMMalgorithmID)

  # LOOK UP ID OF A HARD-CODED GENOTYPE_VERSION
  B73_agpv4_maize282_versionID = find.find_genotype_version(conn, "B73 RefGen_v4_AGPv4_Maize282")
  print("B73 agpv4 maize282 genotype version: ")
  print(B73_agpv4_maize282_versionID)  

  # LOOK UP ID OF A HARD-CODED IMPUTATION_METHOD
  majorAlleleImputationID = find.find_imputation_method(conn, "impute to major allele")
  print("major allele imputation ID: ")
  print(majorAlleleImputationID)  

  # LOOK UP ID OF A HARD-CODED KINSHIP
  kinshipID = find.find_kinship(conn, "/opt/BaxDB/file_storage/kinship_files/4.AstleBalding.synbreed.kinship.csv")
  print("kinshipID: ")
  print(kinshipID)  

  # LOOK UP ID OF A HARD-CODED POPULATION_STRUCTURE
  populationStructureID = find.find_population_structure(conn, "/opt/BaxDB/file_storage/population_structure_files/4.Eigenstrat.population.structure.10PCs.csv")
  print("population structure ID: ")
  print(populationStructureID)

  # PARSE GWAS_RUNS FROM FILE AND ADD TO DB
  insertedGwasRunIDs = insert.insert_gwas_runs_from_gwas_results_file(conn, '/home/mwohl/Downloads/GWASdata/9.mlmmResults.csv', MLMMalgorithmID, B73_agpv4_maize282_versionID, 0.2, 0.2, 0.1, majorAlleleImputationID, kinshipID, populationStructureID)
  print("Inserted gwas_run IDs:")
  print(insertedGwasRunIDs)

  # PARSE GWAS_RESULTS FROM FILE AND ADD TO DB
  insertedGwasResultIDs = insert.insert_gwas_results_from_file(conn, maizeSpeciesID, '/home/mwohl/Downloads/GWASdata/9.mlmmResults.csv', MLMMalgorithmID, 0.2, 0.2, majorAlleleImputationID, B73_agpv4_maize282_versionID, kinshipID, populationStructureID, 0.1)
  print("Inserted gwas result IDs: ")
  print(insertedGwasResultIDs)
