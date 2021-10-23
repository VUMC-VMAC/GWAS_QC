# features and improvements to GWAS QC and TOPMed Imputation Pipeline Documentation
# Vaibhav Janve 2021-10

Proposed/planned improvements/modifications:
* Pre-Imputation QC Part 1
	- add build flag -b 'b38' 
	  allows for usingbuild appropriate 1000G data for PC plots
	- add check for 1000G racefile presence (or code it to be same)
	  currrently the script looks for race file with same prefix as 1000G file stem but doesn't check causing failure much later in script and loss of time.

* Pre-Imputation QC Part 2
	- update the HRC check script to Version 4.3.0  https://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim-v4.3.0.zip
	  NOTE: It is recommended to use v4.3.0 for the TOPMed panel, previous versions will work with the TOPMed panel but will require around 300GB RAM.
	  see https://www.well.ox.ac.uk/~wrayner/tools/ for updates
	- update check perl script
	  add new TOPMed tag or update HRC tag to TOPMed in check perl script to avoid confusion in reference panel used
	
