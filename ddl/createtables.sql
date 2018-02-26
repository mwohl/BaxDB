/* Create the sample table */
DROP TABLE IF EXISTS sample;

CREATE TABLE sample (
    sample_id      SERIAL    PRIMARY KEY,
    sample_name    VARCHAR,
    sample_growout INT,
    sample_line    INT
);


