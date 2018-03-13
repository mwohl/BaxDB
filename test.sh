#! /bin/sh

command -v pg_config > /dev/null 2>&1 || { printf "Command 'pg_config' is required but not found in path. Make sure PostgreSQL client tools are installed. Aborting.\n" 1>&2; exit 1; }

PG_LIBDIR=$(pg_config --pkglibdir)
PG_INSTALLDIR="$PG_LIBDIR/baxdb"

mkdir -p -m 755 "$PG_INSTALLDIR" || { printf "Unable to create directory '$pg_installdir' for installation into PostgreSQL. Aborting.
n" 1>&2 exit 1; }

(
    cd ./c &&
    make &&
    cp array_multi_index.so imputed_genotype.so summarize_variant.so "$PG_INSTALLDIR" &&
    chmod -R 755 "$PG_INSTALLDIR"
)
(
    cd ./lib/tinyint-0.1.1 &&
    make &&
    sed -i -e '1i\\\connect baxdb' -e 's|$libdir\/tinyint|$libdir/baxdb/tinyint|g' tinyint.sql &&
    cp tinyint.so "$PG_INSTALLDIR" &&
    chmod -R 755 "$PG_INSTALLDIR"
)
