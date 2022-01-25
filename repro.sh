#!/bin/sh

# apt-get install libboost-dev
#git clone https://github.com/berthubert/antonie2.git
cd antonie2
git pull
make -j4 skfit gcstats genehisto
mkdir genomes
cd genomes
../../retr.sh auto

echo Starting analysis - will take several hours

tar xzf new_taxdump.tar.gz fullnamelineage.dmp
cd auto

../../gcstats ../fullnamelineage.dmp *.fna.gz
../../skfit
../../genehisto *.fna.gz
cp genomes.csv genomes-$(date +"%Y-%m-%d").csv
# diff againts previous one
# diff -uBb <(sort genomes-2021-10-01.csv | cut -f1,2 -d\;) <(sort genomes-2021-10-11.csv | cut -f1,2 -d\;) | grep ^\+ | cut -f2 -d\;
echo *.fna.gz | xargs -n 5000 -P 4 zgrep ^\> > manifest.txt
cd ../../..
python3 gcskew-article.py
