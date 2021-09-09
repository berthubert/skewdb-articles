#!/bin/sh

# apt-get install libboost-dev
git clone https://github.com/berthubert/antonie2.git
cd antonie2
make -j4 skfit gcstats genehisto
mkdir genomes
cd genomes
tar xf ~/Downloads/genome_assemblies_genome_gff.tar
tar xf ~/Downloads/genome_assemblies_genome_fasta.tar
cd ncbi*
../../gcstats *.fna.gz
../../skfit
../../genehisto *.fna.gz


