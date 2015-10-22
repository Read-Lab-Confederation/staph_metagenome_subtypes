#make file does not convert reports to html at the moment
coverage: ./Data/Final_HMP_Matrix.csv ./staph_metagenome_tools.R
	Rscript -e "library(knitr); knit('./HMP_coverage.Rmd')"
	
covX0.5: ./staph_metagenome_tools.R ./Data/cov0.5
	Rscript -e "library(knitr); knit('./cov0.5X_analysis1.Rmd')"