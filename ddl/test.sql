--CREATE SCHEMA IF NOT EXISTS baxdb;

--\c baxdb

--\connect postgres

--DROP DATABASE baxdb;
--DROP ROLE baxdb_owner;

--CREATE ROLE baxdb_owner WITH
--  LOGIN
--  CREATEROLE
--  ENCRYPTED PASSWORD 'password'
--  ;

--CREATE DATABASE baxdb
--  WITH OWNER = baxdb_owner
--  ENCODING = 'UTF-8'
--  ;

\connect baxdb

CREATE TABLE IF NOT EXISTS species (
  species_id SERIAL PRIMARY KEY,
  shortname VARCHAR(45) UNIQUE NOT NULL,
  binomial VARCHAR(45) NOT NULL,
  subspecies VARCHAR(45),
  variety VARCHAR(45))
  ;
