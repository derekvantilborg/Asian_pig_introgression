## Author: Derek van Tilborg

## Description: Script to perform various analyses on pig SNPs and find introgressed regions from Asian pigs in 
#               European domesticated pigs. This script relies heavily on the SNPRelate package.



# -- prerequisites ------------------------------------------------------------------------------------------------------------

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("SNPRelate")

library(SNPRelate)
library(ggplot2)
library(ggrepel)
library(ggridges)
library(readr)
library(wesanderson)
library(biomaRt)
library(pheatmap)
library(data.table)
library(plyr)
library(stringr)
library(Ckmeans.1d.dp)
library(plotly)

# my go to ggplot theme
my_theme <- theme(legend.position = 'right',
                  text = element_text(size=10),
                  legend.key = element_blank(),
                  legend.background = element_blank(),
                  panel.border = element_blank(),
                  panel.background = element_blank(),
                  plot.title = element_text(hjust = 0.5),
                  panel.grid.major = element_line(colour = "grey95"),
                  panel.grid.minor = element_line(colour = "grey95"))

# -- Funtions ------------------------------------------------------------------------------------------------------------

snp_pca = function(gds_file, meta_data, coverage){
  # It is suggested to use a pruned set of SNPs which are in approximate linkage 
  # equilibrium with each other to avoid the strong influence of SNP clusters in 
  # principal component analysis and relatedness analysis.
  #
  # gds_file: .gds file contianing mutation information
  # meta_data: df containing sample metadata
  # coverage: df with sample.id, coverage
  # returns: dataframe of PC1 and 2 along with sample ids and origins
  
  snps_per_chrom <- unlist(snpgdsLDpruning(gds_file, ld.threshold=0.2))
  snps <- unlist(snps_per_chrom)
  pca <- snpgdsPCA(gds_file, snp.id=snps, num.thread=12)
  perc_var <- pca$varprop*100
  pca_df <- data.frame(sample.id = pca$sample.id,
                       origin = meta_data[[4]][match(pca$sample.id, meta_data[[1]])],
                       species = gsub('_','',meta_data[[2]][match(pca$sample.id, meta_data[[1]])]),
                       coverage = coverage[[2]][match(pca$sample.id, coverage[[1]])],
                       PC1 = pca$eigenvect[,1],
                       PC2 = pca$eigenvect[,2],
                       stringsAsFactors = FALSE)
  return(list(pca_df, perc_var))
}

plot_pca = function(pca_df, pig_type, label = T){
  # Function to plot the pca
  # 
  # pca_df: output from the 'snp_pca' function
  # pig_type: string, name of the pig type
  # label: Bool, plot label

  p = ggplot(pca_df[[1]], aes(x = PC2, y = PC1, label = sample.id, col = species, alpha = coverage)) +
    geom_point(size = 2) +
    ggtitle(paste('PCA plot of biallelic and polymorhpic SNPs in', pig_type)) +
    labs(fill = "Origin", 
         caption = 'Asian Domesticated (n=23), Asian Wild (n=21), European Domesticated (n=20), European Wild (n=21)') +
    scale_x_continuous(breaks = seq(-0.2,0.4,0.1), expand = c(0.05, 0.05))+
    scale_y_continuous(breaks = seq(-0.2,0.4,0.1), expand = c(0.05, 0.05))+
    # scale_x_continuous(breaks=seq(0,1,0.1), limits = c(0,1), expand = c(0,0.01)) +
    # ,limits = c(-0.2,0.4)
    xlab(paste0('PC 2 (', round(pca_df[[2]][2],2), '%)') )+
    ylab(paste0('PC 1 (', round(pca_df[[2]][1],2), '%)') )+ 
    my_theme + 
    theme(panel.border = element_rect( fill = NA),
          legend.position = c(0.85, 0.70))
  if (label){p = p+ geom_text_repel(aes(label=sample.id),hjust=0, vjust=0, show.legend = FALSE)}
  print(p)
}


# ----------- Data --------------------------------------------------------------------------------------------------
## Convert VCF to GDS
# snpgdsVCF2GDS('pigs_all_merged.vcf', "pigs_all_merged.gds", method="biallelic.only")

pigs_merged = snpgdsOpen("pigs_all_merged.gds")
pigs_merged_meta = pigs_all_merged_meta <- read_csv("pigs_all_merged_meta.csv")

# coverage
coverage = read.delim('/Users/derek/Desktop/animal/meta/Depth_coverage_after_mapping.txt', header = F, sep = ' ')
coverage = data.frame(sample.id = coverage$V1, coverage = as.numeric(gsub('X','',coverage$V2)))

# Read GDS file
snp.positions = read.gdsn(index.gdsn(pigs_merged, "snp.position"))
snp.ids = read.gdsn(index.gdsn(pigs_merged, "snp.id"))
snp.chromosomes =  read.gdsn(index.gdsn(pigs_merged, "snp.chromosome"))
snp.alleles =  read.gdsn(index.gdsn(pigs_merged, "snp.allele"))
quals =  read.gdsn(index.gdsn(pigs_merged, "snp.annot/qual"))
samp.id <- read.gdsn(index.gdsn(pigs_merged, "sample.id"))
pop_code = gsub('_',' ',pigs_merged_meta[[2]][match(samp.id, pigs_merged_meta[[1]])])



# ------------ PCA ---------------------------------------------------------------------------------------------

pigs_merged_pca = snp_pca(pigs_merged, pigs_merged_meta, coverage)
plot_pca(pigs_merged_pca, 'all pigs', F)
dev.print(pdf, 'All_pigs_PCA.pdf', height = 5, width = 6)



## ---------------- clustering ---------------------------------------------------------------------------

set.seed(100)
ibs.hc <- snpgdsHCluster(snpgdsIBS(pigs_merged, num.thread=12))
rv <- snpgdsCutTree(ibs.hc, samp.group=as.factor(pop_code))

snpgdsDrawTree(rv, main='Hierarchical clustering of Identity-By-State analysis on SNP genotypes',labels = rv$sample.id, leaflab = "perpendicular")
legend("topright", legend=levels(as.factor(pop_code)), col=c(2,1,3,4), pch=19, ncol=4)
dev.print(pdf, 'All_pigs_IBS_dendogram.pdf', height = 9, width = 16)



## ------ Admixture proportions according to pig ancestries ----------------

# run eigen-analysis
RV <- snpgdsEIGMIX(pigs_merged, num.thread = 12)

# define groups
groups <- list(Asian_dom = samp.id[pop_code == "Asian dom"],
               Asian_wild = samp.id[pop_code == "Asian wild"],
               European_wild = samp.id[pop_code == "European wild"],
               Large_White = samp.id[pop_code == "Large White"])

prop <- snpgdsAdmixProp(RV, groups=groups)
snpgdsAdmixPlot(prop, group=pop_code)



## -------- SNP dissimilarity --------------------------------------------------------

sdm = snpgdsDiss(pigs_merged, num.thread = 12)
row.names(sdm$diss) = gsub('_',' ', sdm$sample.id)
colnames(sdm$diss) = gsub('_',' ', sdm$sample.id)
cat_df = data.frame(group = as.factor(pop_code)); row.names(cat_df) = gsub('_',' ', samp.id)

category        <- c("springgreen", "springgreen4", 'dodgerblue', 'dodgerblue4')
names(category) <- c("Asian dom", "Asian wild", "European wild", "Large White")
anno_colors <- list(group = category)

pheatmap(as.matrix(sdm$diss),
         clustering_method = "complete",
         angle_col = '90',
         border_color = "grey20",
         legend_breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, max(sdm$diss)),
         legend_labels = c('0', '0.2', '0.4', '0.6', '0.8',' 1.0', '1.2', "1 - Bij\n"),
         main = 'Individual dissimilarities based on minor allele frequency and missing rate for each SNP.',
         annotation_col = cat_df,
         annotation_row = cat_df,
         annotation_colors = anno_colors,
         # cluster_rows = F,
         sub = 'asfasfaf')

dev.print(pdf, 'All_pigs_diss_heatmap.pdf', height = 12, width = 14)
# Scale is 1 - βij which is formally described in Weir&Goudet (2017).



## -------- Sample miss rate --------------------------------------------------------

sample_miss_rate = data.frame(sample_id = samp.id, miss_rate = snpgdsSampMissRate(pigs_merged))
sample_miss_rate$origin = pigs_merged_meta[[4]][match(sample_miss_rate$sample_id, pigs_merged_meta[[1]])]
sample_miss_rate$species = gsub('_','',pigs_merged_meta[[2]][match(sample_miss_rate$sample_id, pigs_merged_meta[[1]])])
sample_miss_rate$coverage = coverage[[2]][match(sample_miss_rate$sample_id, coverage[[1]])]
sample_miss_rate$Origin = factor(sample_miss_rate$species, levels = c("Asiandom", "Asianwild", "Europeanwild", "Large White"), ordered = TRUE)
sample_miss_rate$sample_id = factor(sample_miss_rate$sample_id, levels = sample_miss_rate$sample_id, ordered = TRUE)

ggplot(sample_miss_rate, aes(x = miss_rate, y = sample_id, color = Origin, fill = Origin, alpha = log(coverage)))+
  geom_bar(stat = 'identity') +
  coord_flip()+
  labs(title = 'Missing rate of all pig samples.', x = 'Missing rate', y = 'Sample')+
  scale_x_continuous(breaks=seq(0,1,0.1), limits = c(0,1), expand = c(0,0.01)) +
  my_theme + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.print(pdf, 'All_pigs_miss_rate.pdf', height = 6, width = 12)

# check for correlation
cor.test(sample_miss_rate$coverage, sample_miss_rate$miss_rate)



## ---------- Between European Per SNP Fst ------------------------------------------------------------------------

snp_Fst_samples = c(samp.id[pop_code == "Asian dom"], #23
                    samp.id[pop_code == "Asian wild"], #21
                    samp.id[pop_code == "Large White"]) #20
snp_Fst_pop = as.factor(c(rep('Asian', 44), rep('European', 20)))

v <- snpgdsFst(pigs_merged, 
               population= snp_Fst_pop,
               sample.id = snp_Fst_samples, 
               method="W&H02", 
               with.id = T)

# Create a dataframe with chr, snp position and Fst score
Fst_per_SNP = data.frame(chr = as.factor(snp.chromosomes[v$snp.id]),
                         pos = snp.positions[v$snp.id],
                         snp.id = v$snp.id,
                         Fst = v$FstSNP)
Fst_per_SNP$Fst[Fst_per_SNP$Fst < 0] = 0



### ------- Within European Per SNP Fst ------------------------------------

snp_Fst_samples_euro = c(samp.id[pop_code == "European wild"], #21
                          samp.id[pop_code == "Large White"]) #20
snp_Fst_pop_euro = as.factor(c(rep('European Wild', 21), rep('European Domesticated', 20)))

v_euro <- snpgdsFst(pigs_merged, 
               population= snp_Fst_pop_euro,
               sample.id = snp_Fst_samples_euro, 
               method="W&H02", 
               with.id = T)

# Create a dataframe with chr, snp position and Fst score
Fst_per_SNP_euro = data.frame(chr = as.factor(snp.chromosomes[v_euro$snp.id]),
                         pos = snp.positions[v_euro$snp.id],
                         snp.id = v_euro$snp.id,
                         Fst = v_euro$FstSNP)
Fst_per_SNP_euro$Fst[Fst_per_SNP_euro$Fst < 0] = 0



### ----------- Fst density comparison --------------------------------------------------------------------------------------------

# # Only keep positive Fst scores lower than 0.05, discard NAs
Fst_per_SNP_subs = subset(Fst_per_SNP, !is.na(Fst) & Fst < 0.05)
Fst_per_SNP_euro_subs = subset(Fst_per_SNP_euro, !is.na(Fst) & Fst < 0.05)

# Factorise and order chromosomes
Fst_per_SNP_subs$chr <- factor(Fst_per_SNP_subs$chr, levels = rev(as.character(c(1:18))) ,ordered = TRUE)
Fst_per_SNP_euro_subs$chr <- factor(Fst_per_SNP_euro_subs$chr, levels = rev(as.character(c(1:18))) ,ordered = TRUE)

Fst_per_SNP_subs$group = 'Between European'
Fst_per_SNP_euro_subs$group = 'Within European'

# merge both dataframes together
Fst_both = rbind(Fst_per_SNP_euro_subs, Fst_per_SNP_subs)
rm(Fst_per_SNP_subs, Fst_per_SNP_euro_subs)

# density ridge plot of the SNP-wise difference between Fst
ggplot(Fst_both, aes(y=chr)) + 
  geom_density_ridges(aes(x = pos, fill = group , color = group), 
                      alpha = .8) +
  labs(x = "SNP Position",
       y = "Chromosome",
       title = "Position wide density of SNPs with a Fst score < 0.05, comparing genetic variation within and between European pigs",
       caption = "\nBetween population: Asian Domesticated pigs (n = 23) + Asian Wild pigs (n = 21), compared to: Large White European Domesticated (n = 20).\n
      Within population: European Wild pigs (n = 21), compared to: Large White European Domesticated (n = 20).\n
      Fst is calculated by relative Beta estimation (Weir & Hill, 2002), using the implementation from (Buckleton et. al., 2016).") +
  scale_y_discrete(expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(0,max(Fst_per_SNP$pos),10000000), expand = c(0, 0))+#, limits = c(0, max(Fst_per_SNP$pos))) +
  my_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        text = element_text(size=8))

dev.print(pdf, 'FstperSNP_within_between_Eurodom.pdf', height = 6, width = 9)



# ---------- Differential Fst per gene ---------------------------------------------------------------------------------------------

Fst_diff = cbind(Fst_per_SNP, data.frame(Fst_euro = Fst_per_SNP_euro$Fst[match(Fst_per_SNP$snp.id, Fst_per_SNP_euro$snp.id)]))
Fst_diff$Fst_diff = Fst_diff$Fst_euro-Fst_diff$Fst # If the value is high here, the region comes from asian pigs
Fst_diff = subset(Fst_diff, !is.na(Fst_diff))

# GTF file for annotation
gtf <- fread('Desktop/animal/Sus_scrofa.Sscrofa11.1.100.chr.gtf')
setnames(gtf, names(gtf), c("chr","source","type","start","end","score","strand","phase","attributes") )

gtf$gene_id = mgsub::mgsub(grep('gene_id', unlist(str_split(gtf$attributes, ';')), value = T),
             c('gene_id \\"', '\\"'), c("", ""))
gtf$cnt = rep(1, nrow(gtf))

# Summarize GTF file per gene
gtf_summ = ddply(gtf, .(gene_id, chr), summarize,
                 start = min(start),
                 end = max(end),
                 SNP_cnt = sum(cnt))

# Calculate the mean Fst per gene, this can take quite a while
gtf_summ$Fst = rep(NA, nrow(gtf_summ))
for (i in 1:nrow(gtf_summ)){
  cat('\r',(100*i/30477))
  gtf_summ$Fst[i] = mean(  subset(Fst_diff, chr == gtf_summ[i,]$chr & pos >= gtf_summ[i,]$start & pos <= gtf_summ[i,]$end)$Fst_diff  )
}

# remove sex chromosomes (we don't have mutation data on that) and order
gtf_summ = subset(gtf_summ, !chr %in% c('X','Y','MT'))
gtf_summ$chr <- factor(gtf_summ$chr , levels = as.character(1:18) ,ordered = TRUE)

# is there a relation between gene length and Fst? --> No there is not, so there is no need to normalize for gene legth
cor.test((gtf_summ$end - gtf_summ$start), gtf_summ$Fst)

# Histogram
ggplot(gtf_summ, aes(x=Fst)) +
  geom_histogram(fill="white", color="black", bins = 100)+
  geom_vline(aes(xintercept=0), color="dodgerblue2", linetype="dashed")+
  geom_text(x=0.4, y=280, label="Asian genetic influence", hjust = 0) +
  annotate("segment", x = 0.4, y = 250, xend = 0.90, yend = 250, size = 2, arrow = arrow(length = unit(0.1, "inches")),
           colour = "dodgerblue2", lineend = 'butt', linejoin = 'mitre') +
  labs(title="Fst distibution of the difference between genetic variation within and between European pigs", x=expression(Delta~Fst), y = 'Frequency')+
  scale_x_continuous(breaks = seq(-1,1,0.2), expand = c(0, 0))+
  my_theme

dev.print(pdf, 'Fst_diff_hist.pdf', height = 6, width = 10)

# Scatter plot
ggplot(gtf_summ, aes(x = Fst, y = start, col = Fst, alpha = 0.8))+
  geom_point() +
  labs(title='Chromosome-wise mean differential Fst per gene, comparing genetic variation within and between European pigs', 
       x = expression(Delta~Fst), 
       y = 'Start position of gene on chromosome',
       caption = 'A higher differential Fst means high genetic influence of Asian pigs'    )+
  facet_grid(. ~ chr) +
  scale_color_gradientn(colours = wes_palette("Zissou1", 100, type = "continuous"))+
  guides(alpha = FALSE) +
  scale_y_continuous(breaks = seq(0,max(gtf_summ$end),10000000), expand = c(0, 0))+
  my_theme

dev.print(pdf, 'Fst_diff_scatter.pdf', height = 10, width = 30)



# ------- Start making a selection of interesting SNPs using annotation --------------------------------

# subset 2 sigma genes
sig_genes = subset(gtf_summ, Fst >= mean(Fst[!is.na(Fst)]) + 2*sd(Fst[!is.na(Fst)]))

# Find back all SNPs from these genes
# Also take 500bp upstream of genes to capture promotors

Fst_diff_sig = Fst_diff[FALSE,]
Fst_diff_sig$gene_id = ''[F]

for(i in 1:nrow(sig_genes)){
  cat('\r', (100*i)/nrow(sig_genes))
  subs = subset(Fst_diff, chr == as.character(sig_genes[i,]$chr) & pos >= sig_genes[i,]$start-500 & pos <= sig_genes[i,]$end)
  subs$gene_id = sig_genes[i,]$gene_id
  Fst_diff_sig = rbind(Fst_diff_sig, subs) 
}

hist(Fst_diff_sig$Fst_diff, breaks = 100)



# ------- Add VEP scores -------------------------------------------------------

# To find the VEPannotation, all SNPs of interest (Fst_diff_sig) were first saved in VCF format (tsv witg Chrom, pos, ., ref, alt).  
# VEP was then run on this VCF filein offline mode using a Docker container with  
# standard settings. The species and assembly were specifiedas ‘susscrofa’ and ‘Sscrofa11.1’ respectively

VEP_pig <- read_delim("VEP_pig.csv", 
                      "\t", escape_double = FALSE, col_names = FALSE, 
                      col_types = cols_only(X1 = col_guess(), 
                                            X7 = col_guess()), comment = "#", 
                      trim_ws = TRUE)
names(VEP_pig) = c('SNP','VEP')

SNP_names = paste0(snp.chromosomes[Fst_diff_sig$snp.id],'_',
       snp.positions[Fst_diff_sig$snp.id],'_',
       snp.alleles[Fst_diff_sig$snp.id])

Fst_diff_sig$VEP = VEP_pig$VEP[match(SNP_names, VEP_pig$SNP)]
Fst_diff_sig$coding = 'non_coding'
Fst_diff_sig$coding[Fst_diff_sig$VEP %in% c('stop_gained','stop_lost','start_lost','missense_variant','synonymous_variant','stop_retained_variant','coding_sequence_variant')] = 'coding'



# ------- Add pCADD scores -------------------------------------------------------

# pCADDscores are available as supplementary data from Groß et al. and can be found at: http://www.bioinformatics.nl/pCADD/indexed_pPHRED-scores/.

  # Groß, C., Derks, M., Megens, H.-J., Bosse, M., Groenen, M.A.,Reinders, M., De Ridder, D.: 
  # pcadd: Snv prioritisation in sus scrofa.Genetics Selection Evolution52(1), 4 (2020

Fst_diff_sig$pCADD = NA
SNP_names = paste0(snp.chromosomes[Fst_diff_sig$snp.id],'_',
                   snp.positions[Fst_diff_sig$snp.id],'_',
                   gsub('/','_',snp.alleles[Fst_diff_sig$snp.id]))

# I've replaced the first 3 occerences of '\t' to '_' in all pCADD files. This allows for easy matching of SNP names with the very fast match function in R.    
# Next, I've split all pCADD files into files of 10000000 lines because single files were too big.  
# 1)  perl -i.bk -pe '$c = 0; s/\t/_/ while $c++ < 3' -- 1_pCADD-PHRED-scores.tsv
# 2)  split -l 10000000 1_pCADD-PHRED-scores.tsv 1_pCADD-PHRED-scores_

file_list = list.files("/Users/derek/Desktop/animal/pCADD/",pattern="pCADD-PHRED")

for (i in file_list){
  cat('\r',(100*which(file_list == i))/length(file_list))
  pCADD_file <- read_delim(paste0("Desktop/animal/pCADD/",i), 
                           "\t", escape_double = FALSE, col_names = FALSE, 
                           comment = "#", trim_ws = TRUE)
  pCADD_matches = pCADD_file$X2[match(SNP_names, pCADD_file$X1)]
  Fst_diff_sig$pCADD[which(!is.na(pCADD_matches))] = pCADD_matches[!is.na(pCADD_matches)]
}
rm(pCADD_file)

# scatter plot
pmain = ggplot(Fst_diff_sig, aes(y = pCADD, x = Fst_diff, col = coding))+
  geom_point(alpha = 0.25) +
  labs(title = 'Variant-wise differential Fst and pCADD',
       x = expression(Delta~Fst),
       y = 'pCADD score') +
  scale_color_manual(values = c("coding" = "firebrick3",
                                "non_coding" = "dodgerblue3"))+
  my_theme + 
  theme(panel.border = element_rect( fill = NA))

# Marginal densities along x axis
xdens <- axis_canvas(pmain, axis = "x")+
  geom_density(data = Fst_diff_sig, aes(x = Fst_diff, fill = coding),
               alpha = 0.5, size = 0.1) +
  scale_fill_manual(values = c("coding" = "firebrick3",
                                "non_coding" = "dodgerblue3"))

# Marginal densities along y axis
ydens <- axis_canvas(pmain, axis = "y", coord_flip = TRUE)+
  geom_density(data = Fst_diff_sig, aes(x = pCADD, fill = coding),
               alpha = 0.5, size = 0.1)+
  coord_flip() +
  scale_fill_manual(values = c("coding" = "firebrick3",
                               "non_coding" = "dodgerblue3"))

p1 <- insert_xaxis_grob(pmain, xdens, grid::unit(.2, "null"), position = "top")
p2<- insert_yaxis_grob(p1, ydens, grid::unit(.2, "null"), position = "right")
ggdraw(p2)

dev.print(pdf, 'Fst_diff_pCATT_overview.pdf', height = 7, width = 8)



# ------- Find most important genes -------------------------------------------------------

gene_imp = Fst_diff_sig[with(Fst_diff_sig, order(as.numeric(as.character(chr)), pos)),] 

# chromosome start according to SusScrofa11.1
chrom_pos = c(0, 274280751,426136793,558956346,689840209,794364888,965128708,1086972763,1225896210,1365375092,
              1434695013,1513848190,1575438401,1783751528,1925492979,2065875214,2145810227,2209172216)

# Calculate the cumulative chromosome start
gene_imp$cum_pos = NA
for (i in 1:18){gene_imp$cum_pos[gene_imp$chr == i] = chrom_pos[i]}
gene_imp$cum_pos = gene_imp$cum_pos + gene_imp$pos

# K means clustering
gene_imp$cluster = Ckmeans.1d.dp(gene_imp$cum_pos, c(10,100))$cluster

# Scoring the clusters
gene_imp$max_pCADD_per_cluster = NA
for (i in 1:100){gene_imp$max_pCADD_per_cluster[gene_imp$cluster == i] = max(subset(gene_imp, cluster == i)$pCADD)  }
gene_imp$median_diff_Fst_per_cluster = NA
for (i in 1:100){gene_imp$median_diff_Fst_per_cluster[gene_imp$cluster == i] = median(subset(gene_imp, cluster == i)$Fst_diff)  }
gene_imp$combi_score = gene_imp$max_pCADD_per_cluster * gene_imp$median_diff_Fst_per_cluster

# df of all clusters with corresponding combi score
combi_score_per_clust = ddply(gene_imp, .(cluster), summarize, combi_score = mean(combi_score))
combi_score_per_clust = combi_score_per_clust[with(combi_score_per_clust, order(-combi_score)),] 



# ------- Plot the top 10 clusters -------------------------------------------------------

# select the top 10 clusters
top_clusters = subset(gene_imp, cluster %in% combi_score_per_clust[1:10,1])

for (i in unique(top_clusters$gene_id)){
  dropp = which(top_clusters$gene_id == i)
  exclude = which( top_clusters$pCADD[which(top_clusters$gene_id == i)] == max( top_clusters$pCADD[which(top_clusters$gene_id == i)] ) )
  top_clusters$gene_id[dropp[-exclude]] = NA
}

ggplot(top_clusters, aes(x = pos, y = Fst_diff, col = pCADD, shape = coding))+
  geom_hline(yintercept=0,linetype = "dashed", alpha = 0.25) +
  geom_text_repel(aes(label = factor(gene_id))) +
  geom_point(aes(alpha = 0.75)) +
  labs(title='pCADD and differential Fst for high interest variants',
       y = expression(Delta~Fst),
       x = 'Variant position on chromosome',
       caption = 'A higher differential Fst means high genetic influence of Asian pigs.
       High impact variants were selected by clustering all variants, and selecting the 10 clustes that had the highest median differential Fst, while having a mutation with a high pCADD score.'    )+
  scale_color_gradientn(colours = wes_palette("Zissou1", 100, type = "continuous"))+
  guides(alpha = FALSE) +
  facet_grid(chr ~ ., scales = "free_x") + 
  scale_x_continuous(breaks = seq(0,max(top_clusters$pos),5000000), expand = c(0, 0))+
  my_theme +
  theme(panel.border=element_rect(size=1, fill= NA))

dev.print(pdf, 'Fst_diff_scatter_pCADD_high_interest_variants.pdf', height = 15, width = 30)



# ------- Plot the top 10 clusters interactively -------------------------------------------------------

top_clusters = subset(gene_imp, cluster %in% combi_score_per_clust[1:10,1])
top_clusters$allele = snp.alleles[top_clusters$snp.id]

top_clusters = top_clusters[c(1,2,17,3,4,5,6,7,8,9,10,12)]
names(top_clusters)[c(5,6)] = c('Between_Fst', 'Within_Fst')
write.csv(top_clusters, 'Desktop/animal/top_clusters_with_variants_of_interest.csv' , row.names = F)
top_clusters$Chrom = paste0(top_clusters$chr, 
                          '\nPos: ',top_clusters$pos,
                          '\nAllele: ', top_clusters$allele,
                          '\nDelta Fst: ', top_clusters$Fst_diff,
                          '\nGene ID: ', top_clusters$gene_id,
                          '\nVEP: ', top_clusters$VEP,
                          '\nCoding: ', top_clusters$coding,
                          '\npCADD: ', top_clusters$pCADD
                          )

p = ggplot(top_clusters, aes(x = pos, y = Fst_diff, col = pCADD, shape = coding, label = Chrom))+
  geom_hline(yintercept=0,linetype = "dashed", alpha = 0.25) +
  geom_point(aes(alpha = 0.75)) +
  # geom_text(aes(label = allele), alpha = 0) +
  labs(title='pCADD and differential Fst for high interest variants',
       y = 'Delta Fst',
       x = 'Variant position on chromosome',
       caption = 'A higher differential Fst means high genetic influence of Asian pigs.
       High impact variants were selected by clustering all variants, and selecting the 10 clustes that had the highest median differential Fst, while having a mutation with a high pCADD score.'    )+
  scale_color_gradientn(colours = wes_palette("Zissou1", 100, type = "continuous"))+
  guides(alpha = FALSE) +
  facet_grid(chr ~ ., scales = "free_x") + 
  scale_x_continuous(breaks = seq(0,max(top_clusters$pos),5000000), expand = c(0, 0))+
  my_theme +
  theme(panel.border=element_rect(size=1, fill= NA))

ggplotly(
  p = p,
  tooltip=c('Chrom'))








