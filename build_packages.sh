#!/bin/bash

# install_remphasis.sh
# attemps to build all remphasis source packages it can find
# Hanno 2020

cd remphasis && \
rm -f ./src/RcppExports.*
Rscript -e 'library(Rcpp); compileAttributes(".")' && \
cd ..

for i in $(ls -d remphasis*/); do 
    R CMD build ${i%%} 
done
