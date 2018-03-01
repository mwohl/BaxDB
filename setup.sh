#! /bin/sh

#command -v pg_config > /dev/null 2>&1 || { printf "Command 'pg_config' is required but not found in path. Make sure PostgreSQL client tools are installed. Aborting.\n" 1>&2; exit 1; }

#pg_libdir = "$(pg_config --pkglibdir)"
#pg_installdir = "$pg_libdir/BaxDB"

#mkdir -p -m 755 "$pg_installdir" || { printf "Unable to create directory '$pg_installdir' for installation into PostgreSQL. Aborting.\n" 1>&2 exit 1; }
#(
#    cd ./c &&
#    make &&
#    cp array_multi_index.so imputed_genotype.so summarize_variant.so "$pg_installdir" &&
#    chmod -R 755 "$pg_installdir"
#)
#(
#    cd ./lib/tinyint-0.1.1 &&
#    make &&
#    sed -i -e '1i\\\connect BaxDB' -e 's|$libdir\/tinyint|$libdir/BaxDB/tinyint|g' tinyint.sql &&
#    cp tinyint.so "$pg_installdir" &&
#    chmod -R 755 "$pg_installdir"
#)

#sudo -u postgres psql -q -U postgres -f ./ddl/setup.sql || { printf "Unable to perform setup for 'baxdb' database as user 'postgres'. Check UNIX account privileges and pg_hba.conf. Aborting.\n" 1>&2; exit 1; }

sudo -u postgres psql -q -U postgres -f ./ddl/test.sql || { printf "Unable to perform setup for 'baxdb' database as user 'postgres'. Check UNIX account privileges and pg_hba.conf. Aborting.\n" 1>&2; exit 1; }
