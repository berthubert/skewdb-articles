#!/bin/sh

rsync --times --progress -zv rsync://ftp.ncbi.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz . && \
rsync --times -zv rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt bacteria_assembly_summary.txt && \
rsync --times -zv rsync://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt archaea_assembly_summary.txt && \
cat bacteria_assembly_summary.txt archaea_assembly_summary.txt > assembly_summary.txt && \
awk -F "\t" '$12=="Complete Genome" && $11=="latest"{print $20}' assembly_summary.txt > ftpdirpaths && \
awk 'BEGIN{FS=OFS="/";filesuffix="genomic.fna.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' ftpdirpaths | sed s,ftp://,rsync://, >  ftpfilepaths && \
awk 'BEGIN{FS=OFS="/";filesuffix="genomic.gff.gz"}{ftpdir=$0;asm=$10;file=asm"_"filesuffix;print ftpdir,file}' ftpdirpaths | sed s,ftp://,rsync://, >>  ftpfilepaths && \
sort < ftpfilepaths | grep ^rsync:// >  ftpfilepaths.srt && \

cat > wrsync <<EOM
#!/bin/sh
DEST=\$1
shift
rsync --times -v \$@ \$DEST
EOM
chmod +x ./wrsync

cat ftpfilepaths.srt | xargs -P 3 -n 500 ./wrsync $1



