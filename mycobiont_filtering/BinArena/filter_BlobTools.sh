#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=1:mem=10gb
#PBS -J 1-83

#get name of contigs file
ACCESSION=$(head -n $PBS_ARRAY_INDEX /rds/general/user/tbl19/home/macrogen_2022_genomes.txt|tail -n 1)

cd $PBS_O_WORKDIR

#make header for coverage file
echo -e '# contig_id\tread_cov\tbase_cov' > /rds/general/project/theollewellynproject/live/BlobTools/bbmap_alignments/${ACCESSION}_concoct.sam.cov

while read contig; do

 #make filtered coverage file
 grep -w $contig /rds/general/project/theollewellynproject/live/BlobTools/bbmap_alignments/${ACCESSION}.sam.cov >> /rds/general/project/theollewellynproject/live/BlobTools/bbmap_alignments/${ACCESSION}_concoct.sam.cov
 
 #make filtered diamond file
 grep -w $contig /rds/general/project/theollewellynproject/live/BlobTools/diamond_output/${ACCESSION}.vs.uniref90.diamond.taxified.out >> /rds/general/project/theollewellynproject/live/BlobTools/diamond_output/${ACCESSION}_concoct.vs.uniref90.diamond.taxified.out
 
 #make filtered blast file
 grep -w $contig /rds/general/project/theollewellynproject/live/BlobTools/blast_output/${ACCESSION}.vs.Lecanoromycetes_genomes.blastn.taxified.out >>  /rds/general/project/theollewellynproject/live/BlobTools/blast_output/${ACCESSION}_concoct.vs.Lecanoromycetes_genomes.blastn.taxified.out
done < filtered_Asco_contigs/${ACCESSION}_concoct_merged_headers.txt
