#!/bin/bash
#SBATCH --account=def-eadiegwe
#SBATCH --array=1-3
#SBATCH --mem-per-cpu=8G
#SBATCH --time=15:00:00
#SBATCH --mail-user=dayi.li@mail.utoronto.ca
#SBATCH --mail-type=ALL

module load StdEnv/2020 gcc/9.3.0 r/4.2.1 geos gdal udunits
export R_LIBS=~/R/x86_64-pc-linux-gnu-library

R -e 'source("/project/6055726/dli346/GC_count_PP/source/v11acs/v11acs_fit.R")'