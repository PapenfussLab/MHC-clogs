setwd("~/Dropbox/classIs/")

library(ggplot2)
library(plyr)
library(grid)
library(gridExtra)

# Load data
filename = "old/output_old/old/genome_search/opossum/ensembl_R68/custom_hmmer2_genomic.txt"
#filename = "output_genomes/opossum/ensembl_R75/custom_hmmer2_genomic.txt"

column_names = c("domain","chrom","x2","sstart","send","x3","qstart","qend","x4","score","evalue","strand")
dat = read.delim(filename, header=FALSE, col.names=column_names)
dat$position = dat$sstart/1e6
model_names = c("supfam_mhc_fs", "C1-set_fs", "MHC_I_C_fs")
domain_names = c("alpha", "Ig", "C-terminal")
domain_names2 = c("alpha[1]", "alpha[2]", "Ig", "C-terminal")
chroms = c("1","2","3","4","5","6","7","8","X","Un")
dat$domain <- factor(dat$domain, levels=model_names, ordered=TRUE)
dat$chrom = factor(dat$chrom, levels=chroms)
dat$name = as.character(mapvalues(dat$domain, from=model_names, to=domain_names))

index_alpha1 = which(dat$domain=="supfam_mhc_fs" & 0.5*(dat$qstart+dat$qend)<90)
index_alpha2 = which(dat$domain=="supfam_mhc_fs" & 0.5*(dat$qstart+dat$qend)>=90)
dat[index_alpha1,]$name = "alpha[1]"
dat[index_alpha2,]$name = "alpha[2]"

dat$name = factor(dat$name, levels=domain_names2, ordered=TRUE)

p1 = ggplot(subset(dat, score>1), aes(position, score), parse=TRUE) +
      geom_point(size=I(1)) +
      facet_grid(name~chrom, scales="free", space="free_x", labeller=label_parsed) +
      scale_x_continuous(breaks=seq(0,800,by=100)) +
      theme_bw(base_size=11) +
      theme(axis.text.x=element_blank()) +
      theme(axis.ticks=element_blank()) +
      theme(panel.margin=unit(0, "lines")) +
      xlab("Position") + ylab("Score")
print(p1)

x = 38000000
d = subset(dat, score>1 & chrom=="2" & sstart>x-10000 & send<x+10000)
p2 = ggplot(d, aes(sstart, score, label=name), parse=TRUE, labeller=label_parsed) +
    geom_point(parse=TRUE) +
    geom_text(hjust=-0.15) +
    scale_x_continuous(limits=c(x, x+3000)) +
    scale_y_continuous(limits=c(-10, 175)) +
    theme_bw(base_size=11) +
    xlab("Position") + ylab("Score")
print(p2)

pdf("output/Figure_1DE_(draft).pdf", useDingbats=FALSE)
grid.arrange(p1, p2, nrow=2, heights=c(4, 1.5))
dev.off()
