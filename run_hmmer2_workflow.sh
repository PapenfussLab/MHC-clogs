# ---------------------------------------------------------------------------
# Class I prediction pipeline (HMMER3)
#
# Authors: Anthony Papenfuss & Chris Davoren
# 
# Dependencies:
# 	hmmer2
# 	python
# ---------------------------------------------------------------------------

export HMMER2_DIR=/usr/local/bioinfsoftware/hmmer/hmmer-2.3.2
export HMMER3_DIR=/usr/local/bioinfsoftware/hmmer/hmmer3.1-snap20121016.1

export WORK_DIR=$HOME/MHC-clogs
export NUM_CPUs=4
export PATH=$WORK_DIR/bin:$HMMER2_DIR/bin::$HMMER2_DIR/squid:$PATH
export PYTHONPATH=$WORK_DIR/lib:$PYTHONPATH

cd $WORK_DIR

# ---------------------------------------------------------------------------
# Download data
# ---------------------------------------------------------------------------

mhc_setup.py $WORK_DIR/data

# ---------------------------------------------------------------------------
# Perform HMMER2 search using Pfam & custom models
# ---------------------------------------------------------------------------

mkdir $WORK_DIR/output
mkdir $WORK_DIR/output/proteins
mkdir $WORK_DIR/output/proteins/results

search_proteins.py --models hmmer2 $WORK_DIR/data $WORK_DIR/output/proteins
extract_hmmer2_protein_match.py --models hmmer2 $WORK_DIR/data $WORK_DIR/output/proteins
fasta_count $WORK_DIR/output/proteins/results/hmmer2_classI_proteins.fa > $WORK_DIR/output/proteins/results/counts_hmmer2_db_v_custom.txt

hmmalign -q -m data/hmm/hmmer2/MHC_I_fs.hmm \
	$WORK_DIR/output/proteins/results/hmmer2_classI_proteins.fa \
	> data/hmm/custom_hmmer2/MHC_I_fs.sto
hmmalign -q -m data/hmm/hmmer2/MHC_I_ls.hmm \
	$WORK_DIR/output/proteins/results/hmmer2_classI_proteins.fa \
	> data/hmm/custom_hmmer2/MHC_I_ls.sto
hmmalign -q -m data/hmm/hmmer2/supfam_mhc_fs.hmm \
	$WORK_DIR/output/proteins/results/hmmer2_classI_proteins.fa \
	> data/hmm/custom_hmmer2/supfam_mhc_fs.sto
hmmalign -q -m data/hmm/hmmer2/C1-set_fs.hmm \
	$WORK_DIR/output/proteins/results/hmmer2_classI_proteins.fa \
	> data/hmm/custom_hmmer2/C1-set_fs.sto
hmmalign -q -m data/hmm/hmmer2/C1-set_ls.hmm \
	$WORK_DIR/output/proteins/results/hmmer2_classI_proteins.fa \
	> data/hmm/custom_hmmer2/C1-set_ls.sto

hmmbuild -f -F data/hmm/custom_hmmer2/MHC_I_fs.hmm data/hmm/custom_hmmer2/MHC_I_fs.sto
hmmbuild -F data/hmm/custom_hmmer2/MHC_I_ls.hmm data/hmm/custom_hmmer2/MHC_I_ls.sto
hmmbuild -f -F data/hmm/custom_hmmer2/supfam_mhc_fs.hmm data/hmm/custom_hmmer2/supfam_mhc_fs.sto
hmmbuild -f -F data/hmm/custom_hmmer2/C1-set_fs.hmm data/hmm/custom_hmmer2/C1-set_fs.sto
hmmbuild -F data/hmm/custom_hmmer2/C1-set_ls.hmm data/hmm/custom_hmmer2/C1-set_ls.sto

hmmcalibrate data/hmm/custom_hmmer2/MHC_I_fs.hmm
hmmcalibrate data/hmm/custom_hmmer2/MHC_I_ls.hmm
hmmcalibrate data/hmm/custom_hmmer2/supfam_mhc_fs.hmm
hmmcalibrate data/hmm/custom_hmmer2/C1-set_fs.hmm
hmmcalibrate data/hmm/custom_hmmer2/C1-set_ls.hmm

# Search with custom models
search_proteins.py --models custom_hmmer2 $WORK_DIR/data $WORK_DIR/output/proteins
extract_hmmer2_protein_match.py --models custom_hmmer2 $WORK_DIR/data $WORK_DIR/output/proteins
fasta_count $WORK_DIR/output/proteins/results/custom_hmmer2_classI_proteins.fa >> $WORK_DIR/output/proteins/results/counts_hmmer2_db_v_custom.txt
cat $WORK_DIR/output/proteins/results/counts_hmmer2_db_v_custom.txt

#hmmpfam --cpu $NUM_CPUs data/hmm/hmmer2/Pfam_fs $WORK_DIR/output/proteins/results/hmmer2_classI_proteins.fa \
#    > $WORK_DIR/output/proteins/results/hmmer2_classI_proteins_pfam.txt
#hmmpfam --cpu $NUM_CPUs data/hmm/hmmer2/Pfam_fs $WORK_DIR/output/proteins/results/custom_hmmer2_classI_proteins.fa \
#    > $WORK_DIR/output/proteins/results/custom_hmmer2_classI_proteins_pfam.txt

# ---------------------------------------------------------------------------
# Custom hmmer2 genome search
# ---------------------------------------------------------------------------
# Translate 
translate_genomes.py $WORK_DIR/data

# Search genomes
search_genomes.py -t 12 --models custom_hmmer2 $WORK_DIR/data $WORK_DIR/output/genomes > custom_hmmer2.out 2> custom_hmmer2.err &

# Initialise extraction of class I genes
initialize_genome_extract.py -m custom_hmmer2 $WORK_DIR/output/genomes

# Test
extract_hmmer2_genome_match2.py --species opossum --debug -m custom_hmmer2 $WORK_DIR/data $WORK_DIR/output/proteins $WORK_DIR/output/genomes

# Extract class I genes
extract_hmmer2_genome_match2.py --debug -m custom_hmmer2 $WORK_DIR/data $WORK_DIR/output/proteins $WORK_DIR/output/genomes

fasta_rename_merged_seqs.py output/genomes

# ---------------------------------------------------------------------------
