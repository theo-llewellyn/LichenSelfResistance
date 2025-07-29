#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=20:mem=240gb
#PBS -J 1-83

#get name of contigs file
ACCESSION=$(head -n $PBS_ARRAY_INDEX genomes.txt|tail -n 1)
#get the species name
SPECIES=$(head -n $PBS_ARRAY_INDEX filenames.txt|tail -n 1)

module load anaconda3/personal
source activate antismash

cd $PBS_O_WORKDIR

antismash \
 funannotate/${ACCESSION}_CONCOCT_genes_masked/predict_results/${SPECIES}.gbk \
 -c 20 \
 --taxon fungi \
 --cassis \
 --cb-general \
 --cb-subclusters \
 --cb-knownclusters \
 --asf \
 --pfam2go \
 --genefinding-tool none \
 -v \
 --logfile antismash_${ACCESSION}_log.txt \
 --output-dir antismash/${ACCESSION}

conda deactivate
