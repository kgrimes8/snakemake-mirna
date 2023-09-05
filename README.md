# Documentation
### Requirements
Conda
snakemake
Python3

### Pipeline
Snakemake pipeline to extract premirna sequences from fasta and gff3 files, run viennaRNA, and output in tsv format.

Final output is formatted as:
- miRNA ID
- PremiRNA seq (corrected for direction)
- predicted secondary structure in dot bracket notation
- MFE (corrected for temperature)


### To run tests:
Requires conda, access to internet, no previous conda env called 'testenv'
`
cd test-files
bash run-test.py
`