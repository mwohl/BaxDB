-- Connect to the baxdb database
\connect baxdb

-- Create the sample table
DROP TABLE IF EXISTS sample;

--CREATE TABLE IF NOT EXISTS species (
CREATE TABLE species (
  species_id SERIAL PRIMARY KEY,
  shortname VARCHAR(45) UNIQUE NOT NULL,
  binomial VARCHAR(45) NOT NULL,
  subspecies VARCHAR(45),
  variety VARCHAR(45))
  ;
