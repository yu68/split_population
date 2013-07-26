split\_population
================

split\_population based on single nucleosome-level histone mark data
-----

####STEP 1
get absolute intensity of histone modification intensities in single nucleosome level  
_script_: Epi_Intensity_nucleosome.py
```
usage: Epi_Intensity_nucleosome.py [-h] [-N NUCLEOSOME] [-b BAMS [BAMS ...]]
                                   [-f FMT] [-l LEN] [-n NAME [NAME ...]]
                                   [-r RANGES] [-w WEIGHTP [WEIGHTP ...]]
                                   [-o OUTPUT] [-v]

get intensity of epi-modification on single nucleosome

optional arguments:
  -h, --help            show this help message and exit
  -N NUCLEOSOME, --Nucleosome NUCLEOSOME
                        xls file containing nucleosome location information
                        from Danpos output
  -b BAMS [BAMS ...], --bams BAMS [BAMS ...]
                        bed/bam files for epigenetic data
  -f FMT, --fmt FMT     format: bed/bam,default:bam
  -l LEN, --length LEN  average length of ChIP-seq fragment,default:200
  -n NAME [NAME ...], --name NAME [NAME ...]
                        name of each bed sample (to be wrote on the header)
  -r RANGES, --rangeS RANGES
                        search range to find the maximum Epi-intensity in each
                        nucleosome location (default: 100bp)
  -w WEIGHTP [WEIGHTP ...], --weightP WEIGHTP [WEIGHTP ...]
                        parameters to calculate the weight for each
                        read,[half_len of core nucleosome region and half_len
                        of whole regions],default: [75,125]
  -o OUTPUT, --output OUTPUT
                        output file name (can be .txt)
  -v, --verbose         set to output shifted nucleosome centers for each
                        histone mark, one bed file per mark

library dependency: xplib
```

####STEP 2
Generate relative intensities of histone modifications in single nucleosome level from absolute intensities  
_script_: Epi_Intensity_toRelative.R
```
usage: Rscript Epi_Intensity_toRelative.R [-h] [-p PDF] [-o OUTPUT] [-s SIZE]
                                          [-i ITERATION] [--pvalue PVALUE]
                                          input

Generate relative intensities of histone modification in single nucleosome
from absolute intensities

positional arguments:
  input                 input files contain absolute intensities of histone
                        modification in mononucleosome level

optional arguments:
  -h, --help            show this help message and exit
  -p PDF, --pdf PDF     pdf file to store the distribution of maximum
                        intensities for each histone modification. [default
                        "Distribution_of_Max_Intensity.pdf"]
  -o OUTPUT, --output OUTPUT
                        output txt file storing the relative intensities.
                        [default "Epi_Intensity_nucleosome_relative.txt"]
  -s SIZE, --size SIZE  size of the sampled sublist for distribution of max
                        values. [default: 5000]
  -i ITERATION, --iteration ITERATION
                        iteration number to get a pool of max values for its
                        distribution. [default: 1000]
  --pvalue PVALUE       pvalue cutoff of fuzziness and summit value to select
                        high-quality nucleosomes. [default: 0.01]

Need 'argparse' library
``` 
