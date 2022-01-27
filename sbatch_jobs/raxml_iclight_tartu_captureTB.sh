#!/bin/bash -l

# output and error
#SBATCH -o /u/iclight/raxml_phylo_analysis/tartu_tb_capture/slurm.%j.out
#SBATCH -e /u/iclight/raxml_phylo_analysis/tartu_tb_capture/slurm.%j.err
#
# Job Name:
#SBATCH -J raxml_iclight_tartu_captureTB
#
# Directory:
#SBATCH -D /u/iclight/raxml_phylo_analysis/tartu_tb_capture
# Node feature:
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=72
#SBATCH --mem=119000
#SBATCH --mail-type=fail
#SBATCH --mail-type=time_limit
#SBATCH --mail-user=light@mpiib-berlin.mpg.de
#
# wall clock limit:
#SBATCH --time=24:00:00

outdir="/u/iclight/raxml_phylo_analysis/tartu_tb_capture"
conda activate phylo
raxmlHPC-AVX2 -x 060782 -p 060783 -N 250 -f a -T 32 \
-s /ptmp/iclight/eager/mapping/tartu/tb_capture/multi_vcf_analyzer_run/call_conf_20_no_het/results/multivcfanalyzer/snpAlignmentIncludingRefGenome.fasta  \
-m GTRCAT -o canettii.unifiedgenotyper.vcf \
-n tartu_capture_tb_EASI10_tb_phylogeny.nwk -w ${outdir}/
