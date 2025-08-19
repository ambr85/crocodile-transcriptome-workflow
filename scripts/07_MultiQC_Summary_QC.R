#################################
#MULTIQC - RAW DATA (INSTALL SOFTWARE)
#################################
# CREATE a new environment in conda: 
conda create -n multiqc python=3.10
conda activate multiqc

# and then install MultiQC sit the following command:
pip install multiqc==1.14

# Check install with:
multiqc --version
# ➤ Should say: MultiQC v1.14
 
#Select the general file (e.g. rawdata_qc) containing all fastqc reports and then execute:

multiqc .

# This will generate the report of all the files. 

# Reports can also be generated for Mapping STAR log.final reports and for counts after HTSeq. 
