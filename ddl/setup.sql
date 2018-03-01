\connect postgres

DROP DATABASE baxdb;
DROP ROLE baxdb_owner;

CREATE ROLE baxdb_owner WITH
    LOGIN
    CREATEROLE
    ENCRYPTED PASSWORD 'password'
    ;
CREATE DATABASE baxdb
    WITH OWNER = baxdb_owner
    ENCODING = 'UTF-8'
    ;

--\connect baxdb

--CREATE OR REPLACE FUNCTION array_multi_index( ANYARRAY, INTEGER[] )
--    RETURNS ANYARRAY
--    AS '$libdir/baxdb/array_multi_index'
--    LANGUAGE C IMMUTABLE STRICT;
