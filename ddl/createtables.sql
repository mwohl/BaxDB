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
  location_country VARCHAR(45),
  location_state VARCHAR(45),
  location_city VARCHAR(45),
  year INTEGER NOT NULL,
  growout_growout_type INTEGER NOT NULL REFERENCES growout_type (growout_type_id)
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

-- -----------------------
-- Create the sample table
-- -----------------------
DROP TABLE IF EXISTS sample;
CREATE TABLE sample (
  sample_id SERIAL PRIMARY KEY,
  sample_name VARCHAR(45) NOT NULL,
  sample_growout INTEGER NOT NULL REFERENCES growout (growout_id),
  sample_line INTEGER NOT NULL REFERENCES line (line_id)
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

-- ---------------------------
-- Create the nucleotide table
-- ---------------------------
DROP TABLE IF EXISTS nucleotide;
CREATE TABLE nucleotide (
  nucleotide_id SERIAL PRIMARY KEY,
  iupac_nucleotide_code VARCHAR(10) UNIQUE NOT NULL,
  base VARCHAR(25) NOT NULL
  );

-- -------------------------
-- Create the variant table
-- -------------------------
-- untested!
DROP TABLE IF EXISTS variant;
CREATE TABLE variant (
  variant_id SERIAL PRIMARY KEY,
  variant_species INTEGER NOT NULL REFERENCES species (species_id),
  variant_chromosome INTEGER NOT NULL REFERENCES chromosome (chromosome_id),
  variant_pos INTEGER NOT NULL
  );

-- -------------------------
-- Create the genotype table
-- -------------------------
DROP TABLE IF EXISTS genotype;
CREATE TABLE genotype (
  line_ref INTEGER PRIMARY KEY,
  chromosome_id INTEGER NOT NULL,
  genotype tinyint[] NOT NULL,
  FOREIGN KEY (line_ref) REFERENCES line (line_id),
  FOREIGN KEY (chromosome_id) REFERENCES chromosome (chromosome_id)
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
  phenotype_sample INTEGER NOT NULL REFERENCES sample (sample_id),
  phenotype_trait INTEGER NOT NULL REFERENCES trait (trait_id),
  value VARCHAR(45) NOT NULL
  );

-- -----------------------------------
-- Create the reference_genotype table
-- -----------------------------------
DROP TABLE IF EXISTS reference_genotype;
--CREATE TABLE reference_genotype (
--  );

-- -------------------------------
-- Create the gwas_algorithm table
-- -------------------------------
DROP TABLE IF EXISTS gwas_algorithm;
CREATE TABLE gwas_algorithm (
  gwas_algorithm_id SERIAL PRIMARY KEY,
  gwas_algorithm_name VARCHAR(45) UNIQUE NOT NULL
  );

-- ------------------------
-- Create the kinship table
-- ------------------------
DROP TABLE IF EXISTS kinship;
--CREATE TABLE kinship (
--  );

-- ----------------------------------
-- Create the imputation_method table
-- ----------------------------------
DROP TABLE IF EXISTS imputation_method;
CREATE TABLE imputation_method (
  imputation_method_id SERIAL PRIMARY KEY,
  imputation_method_name VARCHAR(45) UNIQUE NOT NULL
  );

-- -------------------------------------
-- Create the population_structure table
-- -------------------------------------
DROP TABLE IF EXISTS population_structure;
--CREATE TABLE population_structure (
--  );

-- -------------------------
-- Create the gwas_run table
-- -------------------------
DROP TABLE IF EXISTS gwas_run;
CREATE TABLE gwas_run (
  gwas_run_id SERIAL PRIMARY KEY,
  gwas_run_gwas_algorithm INTEGER NOT NULL REFERENCES gwas_algorithm (gwas_algorithm_id),
  missing_snp_cutoff_value NUMERIC NOT NULL,
  missing_line_cutoff_value NUMERIC NOT NULL,
  gwas_run_imputation_method INTEGER NOT NULL REFERENCES imputation_method (imputation_method_id)
--  gwas_run_kinship INTEGER NOT NULL REFERENCES kinship (kinship_id),
--  gwas_run_population_structure INTEGER NOT NULL REFERENCES population_structure (population_structure_id)
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
