# ---------------------------------------------------------------------------
# Collect species sequences
# ---------------------------------------------------------------------------

mkdir output
mkdir output/species
for spp_dir in output/genomes/*
do
	spp=$(basename ${spp_dir})
	cp output/genomes/${spp}/ensembl_R75/custom_hmmer2_merged_classI_proteins_renamed.fa output/species/${spp}.fa
done

count_hits_per_species.py > output/counts.txt
excel output/counts.txt

mkdir output/phylo

# ---------------------------------------------------------------------------
# Collect UTs and build trees
# ---------------------------------------------------------------------------

# Manually fix md_UG as this is a short peptide in Ensembl

fasta_extract_UTs.py
clustalo -i output/phylo/UTs.fa -o output/phylo/UTs.aln --outfmt clu --force
clustal2phylip.py output/phylo/UTs.aln output/phylo/UTs.phy

# ../bin/phyml-prottest-macintel -i UTs_edited.phy -d aa -b 500 -m JTT -f e-v e -a e

# ---------------------------------------------------------------------------
# Build tree of selected sequences
# ---------------------------------------------------------------------------

# Compile selected sequences manually

clustalo -i output/selected.fa -o output/selected.aln --outfmt clu --force

# Edit selected.aln in jalview

clustal2phylip.py output/selected_edit.aln output/selected_edit.phy

# Copy to prottest dir

java -jar prottest-3.4.jar -AIC -BIC -F -all-distributions -threads 2 -i examples/MHC_selected_edit.phy -o examples/MHC_selected_edit.txt

# ---------------------------------------------------------------------------
# Build human-opossum trees
# ---------------------------------------------------------------------------

blastp -outfmt 1 -evalue 1e-5 -num_alignments 3 -query data/seq_known/opossum/pep.fa -subject output/species/opossum.fa > output/Ux_vs_opossum.txt
blastp -outfmt 6 -evalue 1e-5 -num_alignments 3 -query data/seq_known/opossum/pep.fa -subject output/species/opossum.fa > output/Ux_vs_opossum6.txt
blastp -outfmt 6 -evalue 1e-5 -num_alignments 3 -subject data/seq_known/opossum/pep.fa -query output/species/opossum.fa > output/Ux_vs_opossum6.txt

cat output/species/human.fa output/species/opossum.fa output/mm_Mills.fa > output/phylo/human_opossum.fa
clustalo -i output/phylo/human_opossum.fa -o output/phylo/human_opossum.aln --outfmt clu --force

# ---------------------------------------------------------------------------
# Build full tree
# ---------------------------------------------------------------------------

cd output/species
cat opossum.fa wallaby.fa devil.fa platypus.fa \
	human.fa mouse.fa dog.fa cow.fa \
	chicken.fa zebrafinch.fa turkey.fa \
	lizard.fa x_tropicalis.fa tetraodon.fa zebrafish.fa \
	c_intestinalis.fa fruitfly.fa lamprey.fa s_cerevisiae.fa > ../phylo/all.fa
cd ../..

clustalo -i output/phylo/all.fa -o output/phylo/all.aln --outfmt clu --force

# EDIT IN JALVIEW & PHYLO WITH BEAST2


sumtrees.py --edges=mean-length --burnin=2000 \
	--support-as-labels --percentages --decimals=1 \
	--min-clade-freq=0.0 --collapse-negative-edges --replace \
	--output=output/phylo/all_edit.nex \
	output/phylo/all1/all_edit.trees \
	output/phylo/all2/all_edit.trees \
	output/phylo/all3/all_edit.trees \
	output/phylo/all4/all_edit.trees


# sumtrees.py --edges=mean-length --burnin=350 \
# 	--support-as-labels --percentages --decimals=1 \
# 	--min-clade-freq=0.15 --replace \
# 	--output=output/all_v1.nex \
# 	output/all_v1/all_edit_nj.trees
#
# sumtrees.py --edges=mean-length --burnin=300 \
# 	--support-as-labels --percentages --decimals=1 \
# 	--min-clade-freq=0.15 --replace \
# 	--output=output/all_v2.nex \
# 	output/all_v2/all_edit_nj.trees
#
#
# sumtrees.py --edges=mean-length --burnin=400 \
# 	--support-as-labels --percentages --decimals=1 \
# 	--min-clade-freq=0.15 --replace \
# 	--output=output/all_v12.nex \
# 	output/all_v1/all_edit_nj.trees \
# 	output/all_v2/all_edit_nj.trees
#
#
# sumtrees.py --edges=mean-length --burnin=900 \
# 	--support-as-labels --percentages --decimals=1 \
# 	--min-clade-freq=0.15 --replace \
# 	--output=output/all_v3_1.nex \
# 	output/all_v3.1/all_edit_nj.trees
#
# sumtrees.py --edges=mean-length --burnin=900 \
# 	--support-as-labels --percentages --decimals=1 \
# 	--min-clade-freq=0.15 --replace \
# 	--output=output/all_v3_2.nex \
# 	output/all_v3.2/all_edit_nj.trees
#
# sumtrees.py --edges=mean-length --burnin=900 \
# 	--support-as-labels --percentages --decimals=1 \
# 	--min-clade-freq=0.15 --replace \
# 	--output=output/all_v3_3.nex \
# 	output/all_v3.3/all_edit_nj.trees
#
# sumtrees.py --edges=mean-length --burnin=800 \
# 	--support-as-labels --percentages --decimals=1 \
# 	--min-clade-freq=0.15 --replace \
# 	--output=output/all_v3.nex \
# 	output/all_v3.1/all_edit_nj.trees \
# 	output/all_v3.3/all_edit_nj.trees
#
# sumtrees.py --edges=mean-length --burnin=0 \
# 	--support-as-labels --percentages --decimals=1 \
# 	--min-clade-freq=0.15 --replace \
# 	--output=output/all_v4.nex \
# 	output/all_v4/all_edit_nj.trees
#
#
# sumtrees.py --edges=mean-length --burnin=16000 \
# 	--support-as-labels --percentages --decimals=1 \
# 	--min-clade-freq=0.15 --replace \
# 	--output=output/all_v5_4000.nex \
# 	output/all_v5/all_edit_nj.trees

# ---------------------------------------------------------------------------
# Compile protein searches
# ---------------------------------------------------------------------------

echo "================================================================================" > output/protein_search.txt
for spp in opossum wallaby devil platypus human mouse dog cow chicken zebrafinch turkey lizard x_tropicalis tetraodon zebrafish c_intestinalis fruitfly lamprey s_cerevisiae
do
	echo $spp >> output/protein_search.txt
	echo "================================================================================" >> output/protein_search.txt
	cat output/proteins/$spp/ensembl_R75/custom_hmmer2_protein_search.txt >> output/protein_search.txt
	echo "================================================================================" >> output/protein_search.txt
done

# ---------------------------------------------------------------------------
# Compile genome searches
# ---------------------------------------------------------------------------

echo "================================================================================" > output/genome_search.txt
for spp in opossum wallaby devil platypus human mouse dog cow chicken zebrafinch turkey lizard x_tropicalis tetraodon zebrafish c_intestinalis fruitfly lamprey s_cerevisiae
do
	echo $spp >> output/genome_search.txt
	echo "================================================================================" >> output/genome_search.txt
	cat output/genomes/$spp/ensembl_R75/custom_hmmer2_chains.txt >> output/genome_search.txt
	echo "================================================================================" >> output/genome_search.txt
done

# ---------------------------------------------------------------------------

