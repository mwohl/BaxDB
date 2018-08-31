# A class is defined per each table in the database

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