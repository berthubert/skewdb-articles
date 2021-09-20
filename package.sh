#!/bin/sh

cp skewdb-readme antonie2/genomes/auto/README
cd antonie2/genomes/auto
tar cjf skewdb.tar.bz2 README gcskewdb.csv skplot.csv *_fit.csv
