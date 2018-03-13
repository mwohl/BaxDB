-- MODIFICATION OF [MySQL Script generated by MySQL Workbench] FOR POSTGRES

---------------------
-- No pgsql equiv? --
---------------------
--SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0;
--SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0;
--SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='TRADITIONAL,ALLOW_INVALID_DATES';

-- -----------------------------------------------------
-- Schema baxdb
-- -----------------------------------------------------

-- -----------------------------------------------------
-- Schema baxdb
-- -----------------------------------------------------
CREATE SCHEMA IF NOT EXISTS baxdb DEFAULT CHARACTER SET utf8 ;
\c baxdb

-- -----------------------------------------------------
-- Table `baxdb`.`species`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS baxdb.species (
  species_id SMALLINT CHECK (species_id > 0) NOT NULL,
  shortname VARCHAR(45) NOT NULL,
  binomial VARCHAR(45) NOT NULL,
  subspecies VARCHAR(45) NULL,
  variety VARCHAR(45) NULL,
  PRIMARY KEY (species_id),
  CONSTRAINT species_id_UNIQUE UNIQUE  (species_id ASC),
  CONSTRAINT shortname_UNIQUE UNIQUE  (shortname ASC))
 ;


-- -----------------------------------------------------
-- Table `baxdb`.`population`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS baxdb.population (
  population_id SMALLINT CHECK (population_id > 0) NOT NULL,
  population_name VARCHAR(100) NOT NULL,
  population_species SMALLINT CHECK (population_species > 0) NOT NULL,
  PRIMARY KEY (population_id),
  CONSTRAINT population_id_UNIQUE UNIQUE  (population_id ASC)
 ,
  CONSTRAINT population_species_species_id
    FOREIGN KEY (population_species)
    REFERENCES baxdb.species (species_id)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
 ;

CREATE INDEX species_id_idx ON baxdb.population (population_species ASC);


-- -----------------------------------------------------
-- Table `baxdb`.`growout_type`
-- -----------------------------------------------------
CREATE SEQUENCE baxdb.growout_type_seq;

CREATE TABLE IF NOT EXISTS baxdb.growout_type (
  growout_type_id SMALLINT CHECK (growout_type_id > 0) NOT NULL DEFAULT NEXTVAL ('baxdb.growout_type_seq'),
  growout_type VARCHAR(45) NOT NULL,
  PRIMARY KEY (growout_type_id),
  CONSTRAINT growout_type_id_UNIQUE UNIQUE  (growout_type_id ASC),
  CONSTRAINT growout_type_UNIQUE UNIQUE  (growout_type ASC))
 ;


-- -----------------------------------------------------
-- Table `baxdb`.`growout`
-- -----------------------------------------------------
CREATE SEQUENCE baxdb.growout_seq;

CREATE TABLE IF NOT EXISTS baxdb.growout (
  growout_id SMALLINT CHECK (growout_id > 0) NOT NULL DEFAULT NEXTVAL ('baxdb.growout_seq'),
  growout_name VARCHAR(75) NOT NULL,
  growout_population SMALLINT CHECK (growout_population > 0) NOT NULL,
  location_country VARCHAR(45) NULL,
  location_state VARCHAR(45) NULL,
  location_city VARCHAR(45) NULL,
  year YEAR NOT NULL,
  growout_growout_type SMALLINT CHECK (growout_growout_type > 0) NOT NULL,
  PRIMARY KEY (growout_id),
  CONSTRAINT growout_id_UNIQUE UNIQUE  (growout_id ASC),
  CONSTRAINT growout_name_UNIQUE UNIQUE  (growout_name ASC)
 ,
  CONSTRAINT growout_growout_type_growout_type_id
    FOREIGN KEY (growout_growout_type)
    REFERENCES baxdb.growout_type (growout_type_id)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT growout_population_population_id
    FOREIGN KEY (growout_population)
    REFERENCES baxdb.population (population_id)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
;

CREATE INDEX growout_growout_type_idx ON baxdb.growout (growout_growout_type ASC);
CREATE INDEX growout_population_idx ON baxdb.growout (growout_population ASC);

-- -----------------------------------------------------
-- Table `baxdb`.`line`
-- -----------------------------------------------------
CREATE SEQUENCE baxdb.line_seq;

CREATE TABLE IF NOT EXISTS baxdb.line (
  line_id SMALLINT CHECK (line_id > 0) NOT NULL DEFAULT NEXTVAL ('baxdb.line_seq'),
  line_name VARCHAR(45) NOT NULL,
  line_population SMALLINT CHECK (line_population > 0) NOT NULL,
  PRIMARY KEY (line_id),
  CONSTRAINT line_id_UNIQUE UNIQUE  (line_id ASC)
 ,
  CONSTRAINT line_name_UNIQUE UNIQUE  (line_name ASC),
  CONSTRAINT line_population_population_id
    FOREIGN KEY (line_population)
    REFERENCES baxdb.population (population_id)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
;

CREATE INDEX population_population_id_idx ON baxdb.line (line_population ASC);


-- -----------------------------------------------------
-- Table `baxdb`.`sample`
-- -----------------------------------------------------
CREATE SEQUENCE baxdb.sample_seq;

CREATE TABLE IF NOT EXISTS baxdb.sample (
  sample_id INT CHECK (sample_id > 0) NOT NULL DEFAULT NEXTVAL ('baxdb.sample_seq'),
  sample_name VARCHAR(45) NOT NULL,
  sample_growout SMALLINT CHECK (sample_growout > 0) NOT NULL,
  sample_line SMALLINT CHECK (sample_line > 0) NOT NULL,
  PRIMARY KEY (sample_id),
  CONSTRAINT sample_id_UNIQUE UNIQUE  (sample_id ASC)
 ,
  CONSTRAINT sample_line_line_id
    FOREIGN KEY (sample_line)
    REFERENCES baxdb.line (line_id)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT sample_growout_growout_id
    FOREIGN KEY (sample_growout)
    REFERENCES baxdb.growout (growout_id)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
;

CREATE INDEX line_line_id_idx ON baxdb.sample (sample_line ASC);
CREATE INDEX growout_growout_id_idx ON baxdb.sample (sample_growout ASC);


-- -----------------------------------------------------
-- Table `baxdb`.`chromosome`
-- -----------------------------------------------------
CREATE SEQUENCE baxdb.chromosome_seq;

CREATE TABLE IF NOT EXISTS baxdb.chromosome (
  chromosome_id SMALLINT CHECK (chromosome_id > 0) NOT NULL DEFAULT NEXTVAL ('baxdb.chromosome_seq'),
  chromosome_species SMALLINT CHECK (chromosome_species > 0) NOT NULL,
  chromosome_name VARCHAR(45) NOT NULL,
  PRIMARY KEY (chromosome_id),
  CONSTRAINT chromosome_id_UNIQUE UNIQUE  (chromosome_id ASC)
 ,
  CONSTRAINT chromosome_species_species_id
    FOREIGN KEY (chromosome_species)
    REFERENCES baxdb.species (species_id)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
;

CREATE INDEX species_species_id_idx ON baxdb.chromosome (chromosome_species ASC);


-- -----------------------------------------------------
-- Table `baxdb`.`nucleotide`
-- -----------------------------------------------------
CREATE SEQUENCE baxdb.nucleotide_seq;

CREATE TABLE IF NOT EXISTS baxdb.nucleotide (
  nucleotide_id SMALLINT CHECK (nucleotide_id > 0) NOT NULL DEFAULT NEXTVAL ('baxdb.nucleotide_seq'),
  iupac_nucleotide_code VARCHAR(10) NOT NULL,
  base VARCHAR(25) NOT NULL,
  PRIMARY KEY (nucleotide_id),
  CONSTRAINT nucleotide_id_UNIQUE UNIQUE  (nucleotide_id ASC),
  CONSTRAINT iupac_nucleotide_code_UNIQUE UNIQUE  (iupac_nucleotide_code ASC))
 ;


-- -----------------------------------------------------
-- Table `baxdb`.`genotype`
-- -----------------------------------------------------
CREATE SEQUENCE baxdb.genotype_seq;

CREATE TABLE IF NOT EXISTS baxdb.genotype (
  genotype_id BIGINT CHECK (genotype_id > 0) NOT NULL DEFAULT NEXTVAL ('baxdb.genotype_seq'),
  genotype_line SMALLINT CHECK (genotype_line > 0) NOT NULL,
  genotype_chromosome SMALLINT CHECK (genotype_chromosome > 0) NOT NULL,
  position INT CHECK (position > 0) NOT NULL,
  genotype_nucleotide SMALLINT CHECK (genotype_nucleotide > 0) NOT NULL,
  PRIMARY KEY (genotype_id),
  CONSTRAINT genotype_id_UNIQUE UNIQUE  (genotype_id ASC)
 ,
  CONSTRAINT genotype_line_line_id
    FOREIGN KEY (genotype_line)
    REFERENCES baxdb.line (line_id)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT genotype_chromosome_chromosome_id
    FOREIGN KEY (genotype_chromosome)
    REFERENCES baxdb.chromosome (chromosome_id)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT genotype_nucleotide_nucleotide_id
    FOREIGN KEY (genotype_nucleotide)
    REFERENCES baxdb.nucleotide (nucleotide_id)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
 ;

CREATE INDEX line_line_id_idx ON baxdb.genotype (genotype_line ASC);
CREATE INDEX chromosome_chromosome_id_idx ON baxdb.genotype (genotype_chromosome ASC);
CREATE INDEX nucleotide_nucleotide_id_idx ON baxdb.genotype (genotype_nucleotide ASC);


-- -----------------------------------------------------
-- Table `baxdb`.`trait`
-- -----------------------------------------------------
CREATE SEQUENCE baxdb.trait_seq;

CREATE TABLE IF NOT EXISTS baxdb.trait (
  trait_id SMALLINT CHECK (trait_id > 0) NOT NULL DEFAULT NEXTVAL ('baxdb.trait_seq'),
  trait_name VARCHAR(45) NOT NULL,
  measurement_unit VARCHAR(45) NULL,
  measurement_device VARCHAR(45) NULL,
  description VARCHAR(255) NULL,
  PRIMARY KEY (trait_id),
  CONSTRAINT trait_id_UNIQUE UNIQUE  (trait_id ASC),
  CONSTRAINT trait_name_UNIQUE UNIQUE  (trait_name ASC))
 ;


-- -----------------------------------------------------
-- Table `baxdb`.`phenotype`
-- -----------------------------------------------------
CREATE SEQUENCE baxdb.phenotype_seq;

CREATE TABLE IF NOT EXISTS baxdb.phenotype (
  phenotype_id INT CHECK (phenotype_id > 0) NOT NULL DEFAULT NEXTVAL ('baxdb.phenotype_seq'),
  phenotype_sample INT CHECK (phenotype_sample > 0) NOT NULL,
  phenotype_trait SMALLINT CHECK (phenotype_trait > 0) NOT NULL,
  value VARCHAR(45) NULL,
  PRIMARY KEY (phenotype_id),
  CONSTRAINT phenotype_id_UNIQUE UNIQUE  (phenotype_id ASC)
 ,
  CONSTRAINT phenotype_sample_sample_id
    FOREIGN KEY (phenotype_sample)
    REFERENCES baxdb.sample (sample_id)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT phenotype_trait_trait_id
    FOREIGN KEY (phenotype_trait)
    REFERENCES baxdb.trait (trait_id)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
 ;

CREATE INDEX sample_sample_id_idx ON baxdb.phenotype (phenotype_sample ASC);
CREATE INDEX trait_trait_id_idx ON baxdb.phenotype (phenotype_trait ASC);


-- -----------------------------------------------------
-- Table `baxdb`.`reference_genotype`
-- -----------------------------------------------------
CREATE TABLE IF NOT EXISTS baxdb.reference_genotype (
  reference_genotype_id INT NOT NULL,
  reference_genotype_name VARCHAR(45) NULL,
  reference_genotype_species VARCHAR(45) NULL,
  reference_genotype_chromosome VARCHAR(45) NULL,
  reference_genotype_position VARCHAR(45) NULL,
  reference_genotype_nucleotide VARCHAR(45) NULL,
  PRIMARY KEY (reference_genotype_id))
 ;


-- -----------------------------------------------------
-- Table `baxdb`.`gwas_algorithm`
-- -----------------------------------------------------
CREATE SEQUENCE baxdb.gwas_algorithm_seq;

CREATE TABLE IF NOT EXISTS baxdb.gwas_algorithm (
  GWAS_algorithm_id SMALLINT CHECK (GWAS_algorithm_id > 0) NOT NULL DEFAULT NEXTVAL ('baxdb.gwas_algorithm_seq'),
  GWAS_algorithm_name VARCHAR(45) NOT NULL,
  PRIMARY KEY (GWAS_algorithm_id),
  CONSTRAINT GWAS_algorithm_id_UNIQUE UNIQUE  (GWAS_algorithm_id ASC),
  CONSTRAINT GWAS_algorithm_name_UNIQUE UNIQUE  (GWAS_algorithm_name ASC))
;


-- -----------------------------------------------------
-- Table `baxdb`.`kinship`
-- -----------------------------------------------------
CREATE SEQUENCE baxdb.kinship_seq;

CREATE TABLE IF NOT EXISTS baxdb.kinship (
  kinship_id SMALLINT CHECK (kinship_id > 0) NOT NULL DEFAULT NEXTVAL ('baxdb.kinship_seq'),
  kinship_calculation_method VARCHAR(45) NULL,
  kinship_matrix_storage_loc VARCHAR(45) NULL,
  PRIMARY KEY (kinship_id),
  CONSTRAINT kinship_id_UNIQUE UNIQUE  (kinship_id ASC))
 ;


-- -----------------------------------------------------
-- Table `baxdb`.`imputation_method`
-- -----------------------------------------------------
CREATE SEQUENCE baxdb.imputation_method_seq;

CREATE TABLE IF NOT EXISTS baxdb.imputation_method (
  imputation_method_id SMALLINT CHECK (imputation_method_id > 0) NOT NULL DEFAULT NEXTVAL ('baxdb.imputation_method_seq'),
  imputation_method_name VARCHAR(45) NOT NULL,
  PRIMARY KEY (imputation_method_id),
  CONSTRAINT imputation_method_id_UNIQUE UNIQUE  (imputation_method_id ASC),
  CONSTRAINT imputation_method_name_UNIQUE UNIQUE  (imputation_method_name ASC))
 ;


-- -----------------------------------------------------
-- Table `baxdb`.`population_structure`
-- -----------------------------------------------------
CREATE SEQUENCE baxdb.population_structure_seq;

CREATE TABLE IF NOT EXISTS baxdb.population_structure (
  population_structure_id SMALLINT CHECK (population_structure_id > 0) NOT NULL DEFAULT NEXTVAL ('baxdb.population_structure_seq'),
  population_structure_calculation_method VARCHAR(45) NULL,
  population_structure_storage_loc VARCHAR(45) NULL,
  PRIMARY KEY (population_structure_id),
  CONSTRAINT population_structure_id_UNIQUE UNIQUE  (population_structure_id ASC))
 ;


-- -----------------------------------------------------
-- Table `baxdb`.`gwas_run`
-- -----------------------------------------------------
CREATE SEQUENCE baxdb.gwas_run_seq;

CREATE TABLE IF NOT EXISTS baxdb.gwas_run (
  GWAS_run_id SMALLINT CHECK (GWAS_run_id > 0) NOT NULL DEFAULT NEXTVAL ('baxdb.gwas_run_seq'),
  gwas_algorithm SMALLINT CHECK (gwas_algorithm > 0) NOT NULL,
  missing_snp_cutoff_value DECIMAL(5,4) NOT NULL,
  missing_line_cutoff_value DECIMAL(5,4) NOT NULL,
  gwas_run_imputation_method SMALLINT CHECK (gwas_run_imputation_method > 0) NOT NULL,
  gwas_run_kinship SMALLINT CHECK (gwas_run_kinship > 0) NOT NULL,
  gwas_run_population_structure SMALLINT CHECK (gwas_run_population_structure > 0) NOT NULL,
  PRIMARY KEY (GWAS_run_id),
  CONSTRAINT GWAS_run_id_UNIQUE UNIQUE  (GWAS_run_id ASC)
 ,
  CONSTRAINT gwas_run_gwas_algorithm_gwas_algorithm_id
    FOREIGN KEY (gwas_algorithm)
    REFERENCES baxdb.gwas_algorithm (GWAS_algorithm_id)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT gwas_run_kinship_kinship_id
    FOREIGN KEY (gwas_run_kinship)
    REFERENCES baxdb.kinship (kinship_id)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT gwas_run_imputation_method_imputation_method_id
    FOREIGN KEY (gwas_run_imputation_method)
    REFERENCES baxdb.imputation_method (imputation_method_id)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT gwas_run_population_structure_population_structure_id
    FOREIGN KEY (gwas_run_population_structure)
    REFERENCES baxdb.population_structure (population_structure_id)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
 ;

CREATE INDEX gwas_algorithm_gwas_algorithm_id_idx ON baxdb.gwas_run (gwas_algorithm ASC);
CREATE INDEX gwas_run_kinship_kinship_id_idx ON baxdb.gwas_run (gwas_run_kinship ASC);
CREATE INDEX imputation_method_imputation_method_id_idx ON baxdb.gwas_run (gwas_run_imputation_method ASC);
CREATE INDEX gwas_run_population_structure_population_structure_id_idx ON baxdb.gwas_run (gwas_run_population_structure ASC);


-- -----------------------------------------------------
-- Table `baxdb`.`gwas_result`
-- -----------------------------------------------------
CREATE SEQUENCE baxdb.gwas_result_seq;

CREATE TABLE IF NOT EXISTS baxdb.gwas_result (
  gwas_result_id INT CHECK (gwas_result_id > 0) NOT NULL DEFAULT NEXTVAL ('baxdb.gwas_result_seq'),
  gwas_result_chromosome SMALLINT CHECK (gwas_result_chromosome > 0) NOT NULL,
  basepair INT CHECK (basepair > 0) NOT NULL,
  gwas_result_trait SMALLINT CHECK (gwas_result_trait > 0) NOT NULL,
  snp DOUBLE PRECISION NULL,
  logp DOUBLE PRECISION NULL,
  null_logp DOUBLE PRECISION NULL,
  cofactor SMALLINT CHECK (cofactor > 0) NULL,
  order SMALLINT CHECK (order > 0) NULL,
  gwas_result_gwas_run SMALLINT CHECK (gwas_result_gwas_run > 0) NOT NULL,
  PRIMARY KEY (gwas_result_id),
  CONSTRAINT gwas_result_id_UNIQUE UNIQUE  (gwas_result_id ASC)
 ,
  CONSTRAINT gwas_result_chromosome_chromosome_id
    FOREIGN KEY (gwas_result_chromosome)
    REFERENCES baxdb.chromosome (chromosome_id)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT gwas_result_trait_trait_id
    FOREIGN KEY (gwas_result_trait)
    REFERENCES baxdb.trait (trait_id)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION,
  CONSTRAINT gwas_result_gwas_run_gwas_run_id
    FOREIGN KEY (gwas_result_gwas_run)
    REFERENCES baxdb.gwas_run (GWAS_run_id)
    ON DELETE NO ACTION
    ON UPDATE NO ACTION)
 ;

CREATE INDEX chromosome_chromosome_id_idx ON baxdb.gwas_result (gwas_result_chromosome ASC);
CREATE INDEX trait_trait_id_idx ON baxdb.gwas_result (gwas_result_trait ASC);
CREATE INDEX gwas_run_gwas_run_id_idx ON baxdb.gwas_result (gwas_result_gwas_run ASC);

---------------------
-- No pgsql equiv? --
---------------------
--SET SQL_MODE=@OLD_SQL_MODE;
--SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS;
--SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS;