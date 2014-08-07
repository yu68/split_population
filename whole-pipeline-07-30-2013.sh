
# Step1: get the absolute intensities of histone marks in single nucleosome level 
#   1. change ChIP-seq fragment length as 500
#   2. alllow search for neiboering regions to find histone mark peaks

folder=~/ChIPseq_map_BAM/

nucleosome_file=/home/yu68/MNase-seq/Danpos/whole_MNase/result/pooled/Pn_E14_whole_mm9.sort.Fnor.ajClonal.smooth.peaks.xls

python Epi_Intensity_nucleosome.py -v -l 300 -N $nucleosome_file -b ${folder}mouse_H3K4me3_d0.sort.bam ${folder}mouse_H3K4me2_d0.sort.bam ${folder}mouse_H3K4me_d0.sort.bam ${folder}mouse_H3K27me3_d0.sort.bam ${folder}mouse_H3K36me3_d0.sort.bam ${folder}mouse_H3K27Ac_d0.sort.bam ${folder}mouse_H2A.Z_d0.sort.bam ${folder}mouse_H3K9me3_d0.sort.bam -n H3K4me3 H3K4me2 H3K4me1 H3K27me3 H3K36me3 H3K27ac H2AZ H3K9me3 -w 75 125 -o Epi_Intensity_nucleosome-300bp_75-125.txt


# Step2: estimate the relative intensities [0-1] from absolute intensities

Rscript Epi_Intensity_toRelative.R -p Distribution_of_Max_Intensity_1E-2_300bp.pdf -o Epi_Intensity_nucleosome_relative_1E-2_300bp.txt Epi_Intensity_nucleosome-300bp_75-125.txt


# Step3: Run the split population program

python split_population_fixedCluster.py Epi_Intensity_nucleosome_relative_1E-2_300bp.txt -c 3 -o _matrix.txt -i 50 -d -v


# Step4: downstream analysis
Rscript plot_colocalization.R Clut-3_Hist-7_Nucleo-974874_beta-500_matrix.txt -n 3 -o colocalization_plot.pdf 


# Step5: visualize examples
python visualize_nucleosome_pattern_chromatinState.py -i Clut-3_Hist-7_Nucleo-974874_beta-500_matrix.txt -r chr6:122651369-122670010 -o visualize_nucleosome_Nanog.pdf


