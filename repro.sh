#!/bin/sh

# apt-get install libboost-dev
git clone https://github.com/berthubert/antonie2.git
cd antonie2
make -j4 skfit gcstats genehisto
mkdir genomes
cd genomes

echo Now download https://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz and press enter
read val


echo Now go to the NCBI download service: 'https://www.ncbi.nlm.nih.gov/assembly/?term=all%5Bfilter%5D'
echo Select '"organism group"' bacteria or archaea, '"assembly level"' 
echo complete genome, and press '"Download Assemblies"'.
echo Then pick '"Genomic FASTA"', and download +- 30GB of tar file
echo After that is done, pick '"Genomic GFF"'  and download another 6.5GB of tar file

echo Once downloaded, put the files in the directory where this script is $(pwd)
read val

echo Unpacking

tar xzf new_taxdump.tar.gz fullnamelineage.dmp
tar xf genome_assemblies_genome_gff.tar
tar xf genome_assemblies_genome_fasta.tar

echo Starting analysis - will take several hours

cd ncbi*
../../gcstats ../../fullnamelineage.dmp *.fna.gz
../../skfit
../../genehisto *.fna.gz


