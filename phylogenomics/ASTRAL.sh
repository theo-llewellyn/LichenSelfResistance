#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=32:mem=124gb

module load anaconda3/personal
source activate iqtree-env 


mkdir ASTRAL/ASTRAL_Leca79T_50p_tAL
cd ASTRAL/ASTRAL_Leca79T_50p_tAL

#make gene trees
iqtree2 -S Leca79T_msa_tAl \
 --prefix gene_OF_tAl \
 -T 32

ASTRAL species tree
#module load openjdk/11 astral/5.7.1
astral -i gene_OF_tAl.treefile \
 -t 3 \
 -o ASTRAL_Leca79T_50p_OF_tAl.tre \
 2>ASTRAL_out.log

conda deactivate
