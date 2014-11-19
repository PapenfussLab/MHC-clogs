# ---------------------------------------------------------------------------
# Class I prediction pipeline (HMMER3)
#
# Authors: Anthony Papenfuss & Chris Davoren
# 
# Dependencies:
#   hmmer3
#   hmmer2
#   python
# ---------------------------------------------------------------------------

export HMMER2_DIR=/usr/local/bioinfsoftware/hmmer/hmmer-2.3.2
export HMMER3_DIR=/usr/local/bioinfsoftware/hmmer/hmmer3.1-snap20121016.1

export WORK_DIR=$HOME/distribution
export NUM_CPUs=4
export PATH=$WORK_DIR/bin:$HMMER3_DIR/bin:$HMMER2_DIR/squid:$PATH
export PYTHONPATH=$WORK_DIR/lib:$PYTHONPATH

cd $WORK

# ---------------------------------------------------------------------------
# Download data
# ---------------------------------------------------------------------------
# mhc_setup.py $WORK/data

# ---------------------------------------------------------------------------
# Perform HMMER3 search using Pfam & custom models
# ---------------------------------------------------------------------------
source $MHC/bashrc_hmmer3

# Turn off acceleration
search_proteins.py --max --models hmmer3 $WORK/data $WORK/output/proteins
extract_hmmer3_protein_match.py --models hmmer3 $WORK/data $WORK/output/proteins
fasta_count $WORK/output/proteins/results/hmmer3_classI_proteins.fa > $WORK/output/proteins/results/counts_hmmer3_db_v_custom.txt

hmmalign --trim $MHC/hmm/hmmer3/supfam_MHC_I.hmm \
	$WORK/output/proteins/results/hmmer3_classI_proteins.fa \
	> $MHC/hmm/custom_hmmer3/supfam_MHC_I.sto
hmmalign --trim $MHC/hmm/hmmer3/MHC_I.hmm \
	$WORK/output/proteins/results/hmmer3_classI_proteins.fa \
	> $MHC/hmm/custom_hmmer3/MHC_I.sto
hmmalign --trim $MHC/hmm/hmmer3/C1-set.hmm \
	$WORK/output/proteins/results/hmmer3_classI_proteins.fa \
	> $MHC/hmm/custom_hmmer3/C1-set.sto

hmmbuild $MHC/hmm/custom_hmmer3/MHC_I.hmm $MHC/hmm/custom_hmmer3/MHC_I.sto
hmmbuild $MHC/hmm/custom_hmmer3/supfam_MHC_I.hmm $MHC/hmm/custom_hmmer3/supfam_MHC_I.sto
hmmbuild $MHC/hmm/custom_hmmer3/C1-set.hmm $MHC/hmm/custom_hmmer3/C1-set.sto

search_proteins.py --max --models custom_hmmer3 $WORK/data $WORK/output/proteins
extract_hmmer3_protein_match.py --models custom_hmmer3 $WORK/data $WORK/output/proteins
fasta_count $WORK/output/proteins/results/custom_hmmer3_classI_proteins.fa >> $WORK/output/proteins/results/counts_hmmer3_db_v_custom.txt
cat $WORK/output/proteins/results/counts_hmmer3_db_v_custom.txt

hmmscan --cpu $NUM_CPUs $MHC/hmm/hmmer3/Pfam-A.hmm $WORK/output/proteins/results/hmmer3_classI_proteins.fa \
    > $WORK/output/proteins/results/hmmer3_classI_proteins_pfam.txt
hmmscan --cpu $NUM_CPUs $MHC/hmm/hmmer3/Pfam-A.hmm $WORK/output/proteins/results/custom_hmmer3_classI_proteins.fa \
    > $WORK/output/proteins/results/custom_hmmer3_classI_proteins_pfam.txt

# ---------------------------------------------------------------------------
# Custom hmmer3 genome search
# ---------------------------------------------------------------------------
translate_genomes.py $WORK/data

search_genomes.py -c -t 8 --max --models hmmer3 $WORK/data $WORK/output/genomes > hmmer3.out 2> hmmer3.err &
search_genomes.py -c -t 8 --max --models custom_hmmer3 $WORK/data $WORK/output/genomes > custom_hmmer3.out 2> custom_hmmer3.err &

initialize_genome_extract.py -m custom_hmmer3 $WORK/output_genome
extract_hmmer3_genome_match.py -m custom_hmmer3 $WORK/data $WORK/output/proteins $WORK/output_genome

cd $WORK/output_genome
cat */ensembl_R68/custom_hmmer3_merged_classI_proteins.fa > phylo/all.fa
cat human/ensembl_R68/custom_hmmer3_merged_classI_proteins.fa \
    mouse/ensembl_R68/custom_hmmer3_merged_classI_proteins.fa \
    opossum/ensembl_R68/custom_hmmer3_merged_classI_proteins.fa \
    wallaby/ensembl_R68/custom_hmmer3_merged_classI_proteins.fa \
    devil/ensembl_R68/custom_hmmer3_merged_classI_proteins.fa \
    platypus/ensembl_R68/custom_hmmer3_merged_classI_proteins.fa \
    > phylo/selected.fa
cat human/ensembl_R68/custom_hmmer3_merged_classI_proteins.fa \
    opossum/ensembl_R68/custom_hmmer3_merged_classI_proteins.fa \
    > phylo/human_opossum.fa

muscle -in human_opossum.fa  -out human_opossum.aln -clwstrict
muscle -in selected.fa  -out selected.aln -clwstrict
muscle -in all.fa  -out all.aln -clwstrict

# ---------------------------------------------------------------------------
