-- Connect to the baxdb database
\connect baxdb

-- ------------------------
-- Create the species table
-- ------------------------
DROP TABLE IF EXISTS species;
-- CREATE TABLE IF NOT EXISTS species (
CREATE TABLE species (
  species_id SERIAL PRIMARY KEY,
  shortname VARCHAR(45) UNIQUE NOT NULL,
  binomial VARCHAR(45) NOT NULL,
  subspecies VARCHAR(45),
  variety VARCHAR(45)
  );

-- ---------------------------
-- Create the population table
-- ---------------------------
DROP TABLE IF EXISTS population;
CREATE TABLE population (
  population_id SERIAL PRIMARY KEY,
  population_name VARCHAR(45) UNIQUE NOT NULL,
  population_species INTEGER REFERENCES species (species_id)
  );

-- -----------------------------
-- Create the growout_type table
-- -----------------------------
DROP TABLE IF EXISTS growout_type;
CREATE TABLE growout_type (
  growout_type_id SERIAL PRIMARY KEY,
  growout_type VARCHAR(45) UNIQUE NOT NULL
  );

-- ------------------------
-- Create the growout table
-- ------------------------
DROP TABLE IF EXISTS growout;
CREATE TABLE growout (
  growout_id SERIAL PRIMARY KEY,
  growout_name VARCHAR(75) UNIQUE NOT NULL,
  growout_population INTEGER NOT NULL REFERENCES population (population_id),
  growout_location INTEGER REFERENCES location (location_id),
  year INTEGER NOT NULL,
  growout_growout_type INTEGER NOT NULL REFERENCES growout_type (growout_type_id)
  );

-- -------------------------
-- Create the location table
-- -------------------------
DROP TABLE IF EXISTS location;
CREATE TABLE location (
  location_id SERIAL PRIMARY KEY,
  country VARCHAR(75) NOT NULL,
  state VARCHAR(75),
  city VARCHAR(75),
  code VARCHAR(2) NOT NULL
  );

-- ---------------------
-- Create the line table
-- ---------------------
DROP TABLE IF EXISTS line;
CREATE TABLE line (
  line_id SERIAL PRIMARY KEY,
  line_name VARCHAR(45) NOT NULL,
  line_population INTEGER NOT NULL REFERENCES population (population_id),
  unique (line_name, line_population)
  );

-- ---------------------------
-- Create the chromosome table
-- ---------------------------
DROP TABLE IF EXISTS chromosome;
CREATE TABLE chromosome (
  chromosome_id SERIAL PRIMARY KEY,
  chromosome_name VARCHAR(45) NOT NULL,
  chromosome_species INTEGER NOT NULL REFERENCES species (species_id),
  unique (chromosome_name, chromosome_species)
  );

-- -------------------------
-- Create the variant table
-- -------------------------
DROP TABLE IF EXISTS variant;
CREATE TABLE variant (
  variant_id SERIAL PRIMARY KEY,
  variant_species INTEGER NOT NULL REFERENCES species (species_id),
  variant_chromosome INTEGER NOT NULL REFERENCES chromosome (chromosome_id),
  variant_pos INTEGER NOT NULL,
  unique (variant_species,variant_chromosome, variant_pos)
  );

-- -------------------------
-- Create the genotype table
-- -------------------------
DROP TABLE IF EXISTS genotype;
CREATE TABLE genotype (
  genotype_id SERIAL PRIMARY KEY,
  genotype_line INTEGER NOT NULL,
  genotype_chromosome INTEGER NOT NULL,
  genotype tinyint[] NOT NULL,
  genotype_genotype_version INTEGER NOT NULL REFERENCES genotype_version (genotype_version_id),
  FOREIGN KEY (genotype_line) REFERENCES line (line_id),
  FOREIGN KEY (genotype_chromosome) REFERENCES chromosome (chromosome_id),
  unique (genotype_line, genotype_chromosome)
  );

-- ----------------------
-- Create the trait table
-- ----------------------
DROP TABLE IF EXISTS trait;
CREATE TABLE trait (
  trait_id SERIAL PRIMARY KEY,
  trait_name VARCHAR(45) UNIQUE NOT NULL,
  measurement_unit VARCHAR(45),
  measurement_device VARCHAR(45),
  description TEXT
  );

-- --------------------------
-- Create the phenotype table
-- --------------------------
DROP TABLE IF EXISTS phenotype;
CREATE TABLE phenotype (
  phenotype_id SERIAL PRIMARY KEY,
  phenotype_line INTEGER NOT NULL REFERENCES line (line_id),
  phenotype_trait INTEGER NOT NULL REFERENCES trait (trait_id),
  phenotype_value VARCHAR(45) NOT NULL
  );

-- -------------------------------
-- Create the gwas_algorithm table
-- -------------------------------
DROP TABLE IF EXISTS gwas_algorithm;
CREATE TABLE gwas_algorithm (
  gwas_algorithm_id SERIAL PRIMARY KEY,
  gwas_algorithm VARCHAR(45) UNIQUE NOT NULL
  );

-- ----------------------------------
-- Create the imputation_method table
-- ----------------------------------
DROP TABLE IF EXISTS imputation_method;
CREATE TABLE imputation_method (
  imputation_method_id SERIAL PRIMARY KEY,
  imputation_method TEXT UNIQUE NOT NULL
  );


-- ----------------------------------
-- Create the kinship_algorithm table
-- ----------------------------------
DROP TABLE IF EXISTS kinship_algorithm;
CREATE TABLE kinship_algorithm (
  kinship_algorithm_id SERIAL PRIMARY KEY,
  kinship_algorithm TEXT NOT NULL
);


-- ------------------------
-- Create the kinship table
-- ------------------------
DROP TABLE IF EXISTS kinship;
CREATE TABLE kinship (
  kinship_id SERIAL PRIMARY KEY,
  kinship_algorithm INTEGER NOT NULL REFERENCES kinship_algorithm (kinship_algorithm_id),
  kinship_file_path TEXT NOT NULL
  );

-- -----------------------------------------------
-- Create the population_structure_algorithm table
-- -----------------------------------------------
DROP TABLE IF EXISTS population_structure_algorith;
CREATE TABLE population_structure_algorithm (
  population_structure_algorithm_id SERIAL PRIMARY KEY,
  population_structure_algorithm TEXT NOT NULL
);


-- -------------------------------------
-- Create the population_structure table
-- -------------------------------------
DROP TABLE IF EXISTS population_structure;
CREATE TABLE population_structure (
  population_structure_id SERIAL PRIMARY KEY,
  population_structure_algorithm INTEGER NOT NULL REFERENCES population_structure_algorithm (population_structure_algorithm_id),
  population_structure_file_path TEXT NOT NULL
  );

-- ---------------------------------
-- Create the genotype_version table
-- ---------------------------------
DROP TABLE IF EXISTS genotype_version;
CREATE TABLE genotype_version (
  genotype_version_id SERIAL PRIMARY KEY,
  genotype_version_name VARCHAR(75) UNIQUE NOT NULL,
  genotype_version VARCHAR(50) NOT NULL,
  reference_genome INTEGER NOT NULL REFERENCES line (line_id),
  genotype_version_population INTEGER NOT NULL REFERENCES population (population_id),
  unique (genotype_version, reference_genome)
  );

-- -------------------------
-- Create the gwas_run table
-- -------------------------
DROP TABLE IF EXISTS gwas_run;
CREATE TABLE gwas_run (
  gwas_run_id SERIAL PRIMARY KEY,
  gwas_run_trait INTEGER NOT NULL REFERENCES trait (trait_id),
  nsnps INTEGER NOT NULL,
  nlines INTEGER NOT NULL,
  gwas_run_gwas_algorithm INTEGER NOT NULL REFERENCES gwas_algorithm (gwas_algorithm_id),
  gwas_run_genotype_version INTEGER NOT NULL REFERENCES genotype_version (genotype_version_id),
  missing_snp_cutoff_value NUMERIC NOT NULL,
  missing_line_cutoff_value NUMERIC NOT NULL,
  minor_allele_frequency_cutoff_value NUMERIC NOT NULL,
  gwas_run_imputation_method INTEGER NOT NULL REFERENCES imputation_method (imputation_method_id)
  gwas_run_kinship INTEGER NOT NULL REFERENCES kinship (kinship_id),
  gwas_run_population_structure INTEGER NOT NULL REFERENCES population_structure (population_structure_id),
  -- not sure if it will work to name a unique constraint in table creation statement like so:
  unique unique_combination_of_all_fields (gwas_run_trait, nsnps, nlines, gwas_run_gwas_algorithm, gwas_run_genotype_version, missing_snp_cutoff_value, missing_line_cutoff_value, minor_allele_frequency_cutoff_value, gwas_run_imputation_method, gwas_run_kinship, gwas_run_population_structure)
  --unique (gwas_run_trait, nsnps, nlines, gwas_run_gwas_algorithm, gwas_run_genotype_version, missing_snp_cutoff_value, missing_line_cutoff_value, minor_allele_frequency_cutoff_value, gwas_run_imputation_method, gwas_run_kinship, gwas_run_population_structure)
  );

-- ----------------------------
-- Create the gwas_result table
-- ----------------------------
DROP TABLE IF EXISTS gwas_result;
CREATE TABLE gwas_result (
  gwas_result_id SERIAL PRIMARY KEY,
  gwas_result_chromosome INTEGER NOT NULL REFERENCES chromosome (chromosome_id),
  basepair INTEGER NOT NULL CHECK (basepair > 0),
  gwas_result_trait INTEGER NOT NULL REFERENCES trait (trait_id),
  snp NUMERIC,
  logp NUMERIC,
  null_logp NUMERIC,
  cofactor NUMERIC CHECK (cofactor > 0),
  _order NUMERIC CHECK (_order > 0),
  gwas_result_gwas_run INTEGER NOT NULL REFERENCES gwas_run (gwas_run_id)
  );
