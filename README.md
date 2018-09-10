# reformatting_sperm_csv_files
Reformatting sperm velocity csv files to make them compatible with sperm_analysis code

sperm_analysis expects separate subfolders within a defined working directory. Each of these subfolders should contain the csv files for the separate populations to be compared (e.g. birds from each different location contained in different subfolders; mice with different haplotypes containd in different subfolders). In each of the csv files, sperm_analysis expects the following column headers to be located in the same column position across files:

Band

VCL

Is.valid

VSL

VAP

(Note, "Band" is equivalent to an individual and/or sample ID).

The code in the current repository reads in a folder full of csv files, creates a working directory (sperm_analysis_working_dir), with subfolders corresponding the variable sperm speed will be compared across (user defined). It will spit out csv files into the appropriate working directories, renaming the columns if necessary. Note, column names should not have spaces in their names.

Usage:
```
reformatting_sperm_csv_files(working_dir,subfolder_variable,Band,VCL,Is.valid,VSL,VAP)
```
Where working_dir = location where csv files are stored

subfolder_variable = the variable in the csv files that the subfolders should be created based upon (e.g. location, haplotype)

Band = what the individual ID/sample ID column is called

VCL = what the VCL column is called

Is.valid = what the Is.valid column is called

VSL = what the VSL column is called

VAP = what the VSL column is called

e.g.
```
reformatting_sperm_csv_files("/Users/alanaalexander/Dropbox/polg_mice/Polg_mice","mtDNA","Sample","VCL","Is valid","VSL","VAP")
```
