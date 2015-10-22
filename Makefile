all: ./Data/Final_HMP_Matrix.csv ./staph_metagenome_tools.R
	Rscript -e "library(knitr); knit('./HMP_coverage.Rmd')"