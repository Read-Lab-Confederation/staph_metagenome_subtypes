#best to use 'make core'

#makefile does not convert reports to html at the moment

all: core 40CC
#40CC
40CC: HMP_coverage40.md cov0.5X_analysis_40.md cov0.025X_analysis_40.md

#core analysis
core: HMP_coverage.md cov0.5X_analysis1.md cov0.025X_analysis1.md unfiltered_analysis1.md NYC_Data.md synthetic_data_plots.md spec_sens_plot.md subtyping_tree_figure.md ST_tree_figure.md

HMP_coverage.md: ./Data/Final_HMP_Matrix.csv ./staph_metagenome_tools.R ./HMP_coverage.Rmd
	Rscript -e "library(knitr); knit('./HMP_coverage.Rmd')"
	
cov0.5X_analysis1.md: ./staph_metagenome_tools.R ./Data/cov0.5 ./cov0.5X_analysis1.Rmd
	Rscript -e "library(knitr); knit('./cov0.5X_analysis1.Rmd')"

cov0.025X_analysis1.md: ./staph_metagenome_tools.R ./Data/cov0.025 ./cov0.025X_analysis1.Rmd
	Rscript -e "library(knitr); knit('./cov0.025X_analysis1.Rmd')"

unfiltered_analysis1.md: ./staph_metagenome_tools.R ./Data/combined ./unfiltered_analysis1.Rmd
	Rscript -e "library(knitr); knit('./unfiltered_analysis1.Rmd')"

NYC_Data.md: ./staph_metagenome_tools.R ./Data/DataTable5-metaphlan-metadata_v19.csv ./Data/runs-to-samples.txt ./Data/FInal_NYC_Staph_MecA_Table.csv ./NYC_Data.Rmd
	Rscript -e "library(knitr); knit('./NYC_Data.Rmd')"
	
synthetic_data_plots.md: ./Data/Microbiome_simulation.txt synthetic_data_plots.Rmd
	Rscript -e "library(knitr); knit('synthetic_data_plots.Rmd')"

spec_sens_plot.md: spec_sens_plot.Rmd
	Rscript -e "library(knitr); knit('spec_sens_plot.Rmd')"

subtyping_tree_figure.md: ~/dm ./staph_metagenome_tools.R ./Data/2114_strain_subtypes.csv subtyping_tree_figure.Rmd ./Data/SNP1752540_sample_tags.txt
	Rscript -e "library(knitr); knit('subtyping_tree_figure.Rmd')"
	
ST398_specific_SNPs.tab ST398-specific-SNPs.md: ~/.staphopia_logon.R ST398-specific-SNPs.Rmd
	Rscript -e "library(knitr); knit('ST398-specific-SNPs.Rmd')"
	
screen-ST398_tongue_dorsum.md: SNPs-in-mpileup.R ./Data/combined ./Data/ST398_specific_SNPs.tab ~/Dropbox/ARTICLES_BY_TDR/2015-staph-metagenome/HMP_PROJECT_SHARED/Mpileup_tongue_dorsum/ screen-ST398_tongue_dorsum.Rmd
	Rscript -e "library(knitr); knit('screen-ST398_tongue_dorsum.Rmd')"
	
HMP_coverage40.md: ./Data/Final_HMP_Matrix.csv ./staph_metagenome_tools.R HMP_coverage40.Rmd
	Rscript -e "library(knitr); knit('./HMP_coverage40.Rmd')"

cov0.5X_analysis_40.md: ./staph_metagenome_tools.R ./Data/cov40_0.5 cov0.5X_analysis_40.Rmd
	Rscript -e "library(knitr); knit('./cov0.5X_analysis_40.Rmd')"
	
cov0.025X_analysis_40.md: ./staph_metagenome_tools.R ./Data/cov40_0.025 cov0.025X_analysis_40.Rmd
	Rscript -e "library(knitr); knit('./cov0.025X_analysis_40.Rmd')"

ST_tree_figure.md: ~/dm ./staph_metagenome_tools.R ./Data/S2_Data.csv ./Data/Final_highest_Beta_real_ST_5.txt ST_tree_figure.Rmd
	Rscript -e "library(knitr); knit('ST_tree_figure.Rmd')"
	
#warning: this script takes a day to run and creates a large cache
staphopia_subtype_phylogeny.md ~/dm: ./Data/strains_used_for_subtype_tests.csv staphopia_subtype_phylogeny.Rmd
	Rscript -e "library(knitr); knit('staphopia_subtype_phylogeny.Rmd')"	
