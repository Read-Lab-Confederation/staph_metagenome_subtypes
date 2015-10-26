#makefile does not convert reports to html at the moment
HMP_coverage.md: ./Data/Final_HMP_Matrix.csv ./staph_metagenome_tools.R
	Rscript -e "library(knitr); knit('./HMP_coverage.Rmd')"
	
cov0.5X_analysis1.md: ./staph_metagenome_tools.R ./Data/cov0.5
	Rscript -e "library(knitr); knit('./cov0.5X_analysis1.Rmd')"

cov0.025X_analysis1.md: ./staph_metagenome_tools.R ./Data/cov0.025
	Rscript -e "library(knitr); knit('./cov0.025X_analysis1.Rmd')"

unfiltered_analysis1.md: ./staph_metagenome_tools.R ./Data/combined
	Rscript -e "library(knitr); knit('./unfiltered_analysis1.Rmd')"

NYC_Data.md: ./staph_metagenome_tools.R ./Data/DataTable5-metaphlan-metadata_v19 ./Data/runs-to-samples.txt ./Data/FInal_NYC_Staph_MecA_Table.csv
	Rscript -e "library(knitr); knit('./NYC_Data.Rmd')"
	
synthetic_data_plots.md: ./Data/Microbiome_simulation.txt
	Rscript -e "library(knitr); knit('synthetic_data_plots.Rmd')"

spec_sens_plot.md: 
	Rscript -e "library(knitr); knit('spec_sens_plot.Rmd')"
	
staphopia_subtype_phylogeny.md ~/dm: ./Data/strains_used_for_subtype_tests.csv
	Rscript -e "library(knitr); knit('staphopia_subtype_phylogeny.Rmd')"	

subtyping_tree_figure.md: ~/dm ./staph_metagenome_tools.R ./Data/2114_strain_subtypes.csv
	Rscript -e "library(knitr); knit('subtyping_tree_figure.Rmd')"