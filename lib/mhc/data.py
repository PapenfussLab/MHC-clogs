#!/usr/bin/env python

"""
mhc.data

Data filenames
"""

import os
import copy


WORK_DIR = os.environ["WORK_DIR"]
HMMER2_DIR = os.environ["HMMER2_DIR"]
HMMER3_DIR = os.environ["HMMER3_DIR"]
NUMBER_HMMER_CPUs = 4
NUMBER_THREADS = 6

hmmer2_model_dir = os.path.join(WORK_DIR, "data/hmm/hmmer2")
hmmer3_model_dir = os.path.join(WORK_DIR, "data/hmm/hmmer3")
custom_hmmer3_model_dir = os.path.join(WORK_DIR, "data/hmm/custom_hmmer3")
custom_hmmer2_model_dir = os.path.join(WORK_DIR, "data/hmm/custom_hmmer2")
exon_hmmer2_model_dir = os.path.join(WORK_DIR, "data/hmm/exons")

# HMMER2 models
hmmer2_models = {
    'MHC_I_fs': os.path.join(hmmer2_model_dir, "MHC_I_fs.hmm"),
    'supfam_mhc_fs': os.path.join(hmmer2_model_dir, "supfam_mhc_fs.hmm"),
    'C1-set_fs': os.path.join(hmmer2_model_dir, "C1-set_fs.hmm"),
    'supfam_ig_fs': os.path.join(hmmer2_model_dir, "supfam_ig_fs.hmm"),
    'MHC_I_C_fs': os.path.join(hmmer2_model_dir, "MHC_I_C_fs.hmm"),
    'MHC_II_alpha_fs': os.path.join(hmmer2_model_dir, "MHC_II_alpha_fs.hmm"),
    'MHC_II_beta_fs': os.path.join(hmmer2_model_dir, "MHC_II_beta_fs.hmm"),
}
hmmer2_protein_models = copy.copy(hmmer2_models)

# Custom HMMER2 models
custom_hmmer2_models = {
    'MHC_I_fs': os.path.join(custom_hmmer2_model_dir, "MHC_I_fs.hmm"),
    'supfam_mhc_fs': os.path.join(custom_hmmer2_model_dir, "supfam_mhc_fs.hmm"),
    'C1-set_fs': os.path.join(custom_hmmer2_model_dir, "C1-set_fs.hmm"),
}
custom_hmmer2_protein_models = copy.copy(custom_hmmer2_models)

# HMMER3 models
hmmer3_models = {
    'MHC_I': os.path.join(hmmer3_model_dir, "MHC_I.hmm"),
    'C1-set': os.path.join(hmmer3_model_dir, "C1-set.hmm"),
    'MHC_I_C': os.path.join(hmmer3_model_dir, "MHC_I_C.hmm"),
    'MHC_II_alpha': os.path.join(hmmer3_model_dir, "MHC_II_alpha.hmm"),
    'MHC_II_beta': os.path.join(hmmer3_model_dir, "MHC_II_beta.hmm"),
}

# Custom HMMER3 models
custom_hmmer3_models = {
    'MHC_I': os.path.join(custom_hmmer3_model_dir, "MHC_I.hmm"),
    'C1-set': os.path.join(custom_hmmer3_model_dir, "C1-set.hmm"),
}

# Exon models
exon_hmmer2_models = {
    'alpha1': os.path.join(exon_hmmer2_model_dir, "alpha1.hmm"),
    'alpha2': os.path.join(exon_hmmer2_model_dir, "alpha2.hmm"),
    'beta1': os.path.join(exon_hmmer2_model_dir, "beta1.hmm"),
    'beta2': os.path.join(exon_hmmer2_model_dir, "beta2.hmm"),
    'C1-set_fs': os.path.join(custom_hmmer2_model_dir, "C1-set_fs.hmm"),
}


# data_species_subdir = "ensembl_R68"
# ensembl_genome_filenames = {
#     'c_intestinalis' : 'release-68/fasta/ciona_intestinalis/dna/Ciona_intestinalis.KH.68.dna.toplevel.fa.gz',
#     'chicken'        : 'release-68/fasta/gallus_gallus/dna/Gallus_gallus.WASHUC2.68.dna.toplevel.fa.gz',
#     'cow'            : 'release-68/fasta/bos_taurus/dna/Bos_taurus.UMD3.1.68.dna.toplevel.fa.gz',
#     'devil'          : 'release-68/fasta/sarcophilus_harrisii/dna/Sarcophilus_harrisii.DEVIL7.0.68.dna.toplevel.fa.gz',
#     'dog'            : 'release-68/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.68.dna.toplevel.fa.gz',
#     'fruitfly'       : 'release-68/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP5.68.dna.toplevel.fa.gz',
#     'human'          : 'release-68/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.68.dna.toplevel.fa.gz',
#     'lamprey'        : 'release-68/fasta/petromyzon_marinus/dna/Petromyzon_marinus.Pmarinus_7.0.68.dna.toplevel.fa.gz',
#     'lizard'         : 'release-68/fasta/anolis_carolinensis/dna/Anolis_carolinensis.AnoCar2.0.68.dna.toplevel.fa.gz',
#     'mouse'          : 'release-68/fasta/mus_musculus/dna/Mus_musculus.GRCm38.68.dna.toplevel.fa.gz',
#     'opossum'        : 'release-68/fasta/monodelphis_domestica/dna/Monodelphis_domestica.BROADO5.68.dna.toplevel.fa.gz',
#     'platypus'       : 'release-68/fasta/ornithorhynchus_anatinus/dna/Ornithorhynchus_anatinus.OANA5.68.dna.toplevel.fa.gz',
#     's_cerevisiae'   : 'release-68/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.EF4.68.dna.toplevel.fa.gz',
#     'tetraodon'      : 'release-68/fasta/tetraodon_nigroviridis/dna/Tetraodon_nigroviridis.TETRAODON8.68.dna.toplevel.fa.gz',
#     'turkey'         : 'release-68/fasta/meleagris_gallopavo/dna/Meleagris_gallopavo.UMD2.68.dna.toplevel.fa.gz',
#     'wallaby'        : 'release-68/fasta/macropus_eugenii/dna/Macropus_eugenii.Meug_1.0.68.dna.toplevel.fa.gz',
#     'x_tropicalis'   : 'release-68/fasta/xenopus_tropicalis/dna/Xenopus_tropicalis.JGI_4.2.68.dna.toplevel.fa.gz',
#     'zebrafinch'     : 'release-68/fasta/taeniopygia_guttata/dna/Taeniopygia_guttata.taeGut3.2.4.68.dna.toplevel.fa.gz',
#     'zebrafish'      : 'release-68/fasta/danio_rerio/dna/Danio_rerio.Zv9.68.dna.toplevel.fa.gz',
# }

# data_species_subdir = "ensembl_R74"
# ensembl_genome_filenames = {
#     'c_intestinalis' : 'release-74/fasta/ciona_intestinalis/dna/Ciona_intestinalis.KH.74.dna.toplevel.fa.gz',
#     'chicken'        : 'release-74/fasta/gallus_gallus/dna/Gallus_gallus.Galgal4.74.dna.toplevel.fa.gz',
#     'cow'            : 'release-74/fasta/bos_taurus/dna/Bos_taurus.UMD3.1.74.dna.toplevel.fa.gz',
#     'devil'          : 'release-74/fasta/sarcophilus_harrisii/dna/Sarcophilus_harrisii.DEVIL7.0.74.dna.toplevel.fa.gz',
#     'dog'            : 'release-74/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.74.dna.toplevel.fa.gz',
#     'fruitfly'       : 'release-74/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP5.74.dna.toplevel.fa.gz',
#     'human'          : 'release-74/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.74.dna.toplevel.fa.gz',
#     'lamprey'        : 'release-74/fasta/petromyzon_marinus/dna/Petromyzon_marinus.Pmarinus_7.0.74.dna.toplevel.fa.gz',
#     'lizard'         : 'release-74/fasta/anolis_carolinensis/dna/Anolis_carolinensis.AnoCar2.0.74.dna.toplevel.fa.gz',
#     'mouse'          : 'release-74/fasta/mus_musculus/dna/Mus_musculus.GRCm38.74.dna.toplevel.fa.gz',
#     'opossum'        : 'release-74/fasta/monodelphis_domestica/dna/Monodelphis_domestica.BROADO5.74.dna.toplevel.fa.gz',
#     'platypus'       : 'release-74/fasta/ornithorhynchus_anatinus/dna/Ornithorhynchus_anatinus.OANA5.74.dna.toplevel.fa.gz',
#     's_cerevisiae'   : 'release-74/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.EF4.74.dna.toplevel.fa.gz',
#     'tetraodon'      : 'release-74/fasta/tetraodon_nigroviridis/dna/Tetraodon_nigroviridis.TETRAODON8.74.dna.toplevel.fa.gz',
#     'turkey'         : 'release-74/fasta/meleagris_gallopavo/dna/Meleagris_gallopavo.UMD2.74.dna.toplevel.fa.gz',
#     'wallaby'        : 'release-74/fasta/macropus_eugenii/dna/Macropus_eugenii.Meug_1.0.74.dna.toplevel.fa.gz',
#     'x_tropicalis'   : 'release-74/fasta/xenopus_tropicalis/dna/Xenopus_tropicalis.JGI_4.2.74.dna.toplevel.fa.gz',
#     'zebrafinch'     : 'release-74/fasta/taeniopygia_guttata/dna/Taeniopygia_guttata.taeGut3.2.4.74.dna.toplevel.fa.gz',
#     'zebrafish'      : 'release-74/fasta/danio_rerio/dna/Danio_rerio.Zv9.74.dna.toplevel.fa.gz',
# }


data_species_subdir = "ensembl_R75"
ensembl_genome_filenames = {
    'c_intestinalis' : 'release-75/fasta/ciona_intestinalis/dna/Ciona_intestinalis.KH.75.dna.toplevel.fa.gz',
    'chicken'        : 'release-75/fasta/gallus_gallus/dna/Gallus_gallus.Galgal4.75.dna.toplevel.fa.gz',
    'cow'            : 'release-75/fasta/bos_taurus/dna/Bos_taurus.UMD3.1.75.dna.toplevel.fa.gz',
    'devil'          : 'release-75/fasta/sarcophilus_harrisii/dna/Sarcophilus_harrisii.DEVIL7.0.75.dna.toplevel.fa.gz',
    'dog'            : 'release-75/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.75.dna.toplevel.fa.gz',
    'fruitfly'       : 'release-75/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP5.75.dna.toplevel.fa.gz',
    'human'          : 'release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.toplevel.fa.gz',
    'lamprey'        : 'release-75/fasta/petromyzon_marinus/dna/Petromyzon_marinus.Pmarinus_7.0.75.dna.toplevel.fa.gz',
    'lizard'         : 'release-75/fasta/anolis_carolinensis/dna/Anolis_carolinensis.AnoCar2.0.75.dna.toplevel.fa.gz',
    'mouse'          : 'release-75/fasta/mus_musculus/dna/Mus_musculus.GRCm38.75.dna.toplevel.fa.gz',
    'opossum'        : 'release-75/fasta/monodelphis_domestica/dna/Monodelphis_domestica.BROADO5.75.dna.toplevel.fa.gz',
    'platypus'       : 'release-75/fasta/ornithorhynchus_anatinus/dna/Ornithorhynchus_anatinus.OANA5.75.dna.toplevel.fa.gz',
    's_cerevisiae'   : 'release-75/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.75.dna.toplevel.fa.gz',
    'tetraodon'      : 'release-75/fasta/tetraodon_nigroviridis/dna/Tetraodon_nigroviridis.TETRAODON8.75.dna.toplevel.fa.gz',
    'turkey'         : 'release-75/fasta/meleagris_gallopavo/dna/Meleagris_gallopavo.UMD2.75.dna.toplevel.fa.gz',
    'wallaby'        : 'release-75/fasta/macropus_eugenii/dna/Macropus_eugenii.Meug_1.0.75.dna.toplevel.fa.gz',
    'x_tropicalis'   : 'release-75/fasta/xenopus_tropicalis/dna/Xenopus_tropicalis.JGI_4.2.75.dna.toplevel.fa.gz',
    'zebrafinch'     : 'release-75/fasta/taeniopygia_guttata/dna/Taeniopygia_guttata.taeGut3.2.4.75.dna.toplevel.fa.gz',
    'zebrafish'      : 'release-75/fasta/danio_rerio/dna/Danio_rerio.Zv9.75.dna.toplevel.fa.gz',
}


species_short = {
    'c_intestinalis' : 'ci',
    'chicken'        : 'gg',
    'cow'            : 'bt',
    'devil'          : 'sh',
    'dog'            : 'cf',
    'fruitfly'       : 'dm',
    'human'          : 'hs',
    'lamprey'        : 'pm',
    'lizard'         : 'ac',
    'mouse'          : 'mm',
    'opossum'        : 'md',
    'platypus'       : 'oa',
    's_cerevisiae'   : 'sc',
    'tetraodon'      : 'tn',
    'turkey'         : 'mg',
    'wallaby'        : 'me',
    'x_tropicalis'   : 'xt',
    'zebrafinch'     : 'tg',
    'zebrafish'      : 'dr'
}


ensembl_protein_filenames = {}
for species in sorted(ensembl_genome_filenames.keys()):
    ensembl_protein_filename = ensembl_genome_filenames[species].replace("dna", "pep")
    ensembl_protein_filename = ensembl_protein_filename.replace("toplevel", "all")
    ensembl_protein_filenames[species] = ensembl_protein_filename

ensembl_transcript_filenames = {}
for species in sorted(ensembl_genome_filenames.keys()):
    ensembl_transcript_filename = ensembl_genome_filenames[species].replace("dna", "cdna")
    ensembl_transcript_filename = ensembl_transcript_filename.replace("toplevel", "all")
    ensembl_transcript_filenames[species] = ensembl_transcript_filename


protein_filenames = {}
for species in sorted(ensembl_genome_filenames.keys()):
    basename = os.path.splitext(os.path.basename(ensembl_protein_filenames[species]))[0]
    protein_filenames[species] = os.path.join(species, data_species_subdir, basename)

transcript_filenames = {}
for species in sorted(ensembl_genome_filenames.keys()):
    basename = os.path.splitext(os.path.basename(ensembl_transcript_filenames[species]))[0]
    transcript_filenames[species] = os.path.join(species, data_species_subdir, basename)

genome_filenames = {}
for species in sorted(ensembl_genome_filenames.keys()):
    basename = os.path.splitext(os.path.basename(ensembl_genome_filenames[species]))[0]
    genome_filenames[species] = os.path.join(species, data_species_subdir, basename)

biomart_filenames = {}
for species in sorted(ensembl_genome_filenames.keys()):
    basename = os.path.splitext(os.path.basename(ensembl_genome_filenames[species]))[0]
    biomart_filenames[species] = os.path.join(species, data_species_subdir, 'biomart_export_tsv.txt')
    