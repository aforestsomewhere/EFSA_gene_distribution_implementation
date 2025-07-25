Implementation of "Annex D – Pipeline for the automated analysis of gene distribution in microbial species", available at https://zenodo.org/records/12608405, for HPC environments where Docker/privileges are not available, along with minor bug fixes.

Step 1: Create a conda environment containing the required software
```
mamba create -n efsa_distro_env -c bioconda pip
source activate efsa_distro_env
pip install biopython aniclustermap requests appdirs tqdm
```
Step 2: Load necessary environment + software (modules loaded in this implementation could be rolled into the conda environment if desired)
```
source activate efsa_distro_env
module load fastani/1.33
module load blast/2.14.1
```
Step 3: Run pipeline
Sample command:
```
#Species where sufficient number of complete genomes is available
./run_pipeline.sh test_paralvei "Hafnia paralvei" alveicin_gene.fasta nucleotide extract complete
#Species where few complete genomes are available
./run_pipeline.sh test_paralvei "Hafnia paralvei" alveicin_gene.fasta nucleotide remove incomplete
```
Step 4: Modify cutoffs (if necessary)
Default cutoffs are set at 70% ID and 70% coverage. These can be amended in the main_script.sh script (python3 graphical_out.py)

Note: July 2025 added option to include incomplete genomes. To enable visualisation, have added auxiliary script concatenate_fasta.py, which converts incomplete multifasta genomic data to a single fasta file for each strain.
