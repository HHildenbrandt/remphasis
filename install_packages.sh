#~/bin/bash

# install_remphasis.sh
# attemps to installs remphasis source packages it can find
# Hanno 2020

for i in $(ls remphasis*.gz); do
    Rscript -e "install.packages('${i}', repos = NULL)"
done
