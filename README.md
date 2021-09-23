# Chen-Velikov
 Code used to create results in Chen and Velikov (WP, 2021), Zeroing in on the expected returns of anomalies 

This repository contains code used to create the results in Chen and Velikov (WP, 2021), Zeroing in on the expected returns of anomalies. This code is to be used in conjunction with the MATLAB asset pricing package that accompanies Novy-Marx and Velikov (WP, 2021), Assaying Anomalies. 

The order of operations to replicate the results in Chen and Velikov (WP, 2021) is:

1. Download the code and follow the instructions for setting up the MATLAB asset pricing package from https://github.com/velikov-mihail/AnomalyCookbookOfficial
 	* The results in Chen and Velikov (2021) use the pre-release v0.1 of the MATLAB asset pricing package, which will made public soon. The code is also available from the authors upon request.
2. Download the following two files from Andrew Chen and Tom Zimmerman's www.openassetpricing.com (April 2021 release)
	* [signed_predictors_dl_wide.zip](https://drive.google.com/file/d/1-1RUq2wUADu_ncvQJCYY3wxDhqhHxBow/view) - unzip to get signed_predictors_dl_wide.csv
	* [SignalDocumentation.xlsx](https://docs.google.com/spreadsheets/d/18DvZPscKsD0_ZeeUMjyXhF1qn0emDVaj/edit#gid=70837236)
3. Download the code and follow the instructions for calculating the high-frequency effective spreads from TAQ and ISSM from https://github.com/chenandrewy/hf-spreads-all
 	* After running the code on the WRDS cloud, you should download the output file, hf_monthly.csv
5. Download the code in this repository.
6. Open chen_velikov_main.m and run each cell at a time. The script requires setting up the directories for the MATLAB asset pricing package repository, this repository, and a folder containing the three input files (signed_predictors_dl_wide.csv, SignalDocumentation.xlsx, hf_monthly.csv). It calls multiple other scripts which perform the following functions:  
	* make_tcost_measures.m creates a structure with the individual trading cost measures as well as the composite measure plotted in Figure 3
	* run_unmitigated_strategies.m creates and stores a structure with portfolio sort results without any trading cost mitigation
	* run_mitigated_strategies.m creates and stores several structures with portfolios sort results with trading cost mitigation techniques
	* run_nmv_reconciliation_results.m creates and stores several structures with portfolio sort results reported in Table 6
	* run_combination_strategies.m creates and stores expected return matrices that are constructed by combining predictive signals, as well as a structure with portfolio sort results reported in Table 5
	* organize_results.m organizes the results and stores several structures that combine the unmitigated and mitigated portfolio sorts 
	* make_tables.m prints latex output for all tables in the paper and stores it in Results/ChenVelikovTablesOutput.txt
	* make_figures.m creates all figures in the paper and stores them as .pdf files in /Figures/
	* replicate_tcost_figures.m replicates figures in the original papers for three of the individual trading cost measures (Gibbs, HL, CHL).  
   
