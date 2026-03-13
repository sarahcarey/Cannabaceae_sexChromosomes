###### R script in support of Carey et al., 2026 "An X-linked sex determination mechanism in cannabis and hop"


################# MAIN FIGURES #################

##### Figure 1a, sex chromosome Y-mers #####

library(karyoploteR)
#v1.32.0

## read in data
sex_chroms <- read.table("figure_data/Cannabaceae.sexChroms.lengths.karyoploter.txt", header=TRUE)
sex_chroms.GR <- toGRanges(data.frame(chr=sex_chroms$Chromosome, start=sex_chroms$ChromStart, 
                                  end=sex_chroms$ChromEnd))

coverage_data_Ymers_hops <- read.csv("figure_data/hops_Ymers_1MB.csv", header=TRUE)
coverage_data_Ymers_hops.GR <- toGRanges(data.frame(chr=coverage_data_Ymers_hops$Chromosome, 
                                               start=coverage_data_Ymers_hops$Start, end=coverage_data_Ymers_hops$Stop))

coverage_data_Ymers_cannabis <- read.csv("figure_data/cannabis_Ymers_1MB.csv", header=TRUE)
coverage_data_Ymers_cannabis.GR <- toGRanges(data.frame(chr=coverage_data_Ymers_cannabis$Chromosome, 
                                                    start=coverage_data_Ymers_cannabis$Start, end=coverage_data_Ymers_cannabis$Stop))


## adjust plot parameters

pp <- getDefaultPlotParams(plot.type = 1)
pp$ideogramheight <- 100
pp$data2height <- 0
pp$data1height <- 50
pp$data1outmargin <- 0
pp$data2outmargin <- 0


## generate plot

pdf("figure_output/YmerCov_SexChroms_plot.pdf")

kp <- plotKaryotype(sex_chroms.GR, pin=8, plot.type = 1, labels.plotter = NULL, plot.params = pp)
kpAddChromosomeNames(kp, xoffset=0, cex=0.75)
kpAddBaseNumbers(kp, tick.dist = 25000000, tick.len = 10, tick.col="black", cex=0.75,
                 minor.tick.dist = 10000000, minor.tick.len = 10, minor.tick.col = "black", add.units = F)

kp <-kpHeatmap(kp, data=coverage_data_Ymers_hops.GR, data.panel = "ideogram", 
               y=coverage_data_Ymers_hops$Ymers, colors=c("black", "seagreen1"))

kp <-kpHeatmap(kp, data=coverage_data_Ymers_cannabis.GR, data.panel = "ideogram", 
               y=coverage_data_Ymers_cannabis$Ymers, colors=c("black", "seagreen1"))

dev.off()

##### Figure 1b, GENESPACE gene/repeat/synteny #####
## analysis by John Lovell


##### Figure 1c, GENESPACE gene/repeat/synteny for Cannabis Otto II #####


library("GENESPACE")
#v1.3.1

filepaths <- list(
  assemblies = c(
    "/Users/sarahcarey/Desktop/hemp_genomes/final_assemblies/Cannabis_sativa_OttoII_HAP2_v1.0.X.fa",
    "/Users/sarahcarey/Desktop/hemp_genomes/final_assemblies/Cannabis_sativa_OttoII_HAP1_v1.0.Y.fa"),
  genes = 
    c("/Users/sarahcarey/Desktop/hemp_genomes/annotations/OIIb.primary_high_confidence_ChrRenamed.X.gff3",
    "/Users/sarahcarey/Desktop/hemp_genomes/annotations/OIIa.primary_high_confidence_ChrRenamed.Y.gff3"),
  reps = 
    c("/Users/sarahcarey/Desktop/hemp_genomes/annotations/repeats/EDTA/Cannabis_sativa_OttoII_HAP2_v1.0.fa.mod.EDTA.TEanno.X.gff3",
    "/Users/sarahcarey/Desktop/hemp_genomes/annotations/repeats/EDTA/Cannabis_sativa_OttoII_HAP1_v1.0.fa.mod.EDTA.TEanno.Y.gff3"))

Ymers <- Biostrings::readDNAStringSet("/Users/sarahcarey/Desktop/hemp_genomes/final_assemblies/Cannabis_only_male_kmers.txt.fasta")
Ymers <- as.character(Ymers)

test <- plot_2genomes(
  wd = "/Users/sarahcarey/Desktop/hemp_genomes/genespace",
  genomeIDs <- c("Otto_X", "Otto_Y"),
  faFiles = filepaths$assemblies,
  geneGffFiles = filepaths$genes,
  repeatGffFiles = filepaths$reps,
  nCores = 12, repeatClassColumnName = "type",
  repeatGrep1 = "Gypsy",
  repeatGrep2 = "Copia", returnSourceData = TRUE,
  overwrite = TRUE,
  forceCleanWindows = FALSE, kmers = Ymers)



##### Figure 1d, GENESPACE gene/repeat/synteny for Humulus 21110M #####


library("GENESPACE")
#v1.3.1

filepaths <- list(
  assemblies = c(
    "/Users/sarahcarey/Desktop/hop_genomes/final_assemblies/Humulus_lupulus_lupulus_21110M_HAP2_v1.0.ChrX.fa",
    "/Users/sarahcarey/Desktop/hop_genomes/final_assemblies/Humulus_lupulus_lupulus_21110M_HAP1_v1.2.ChrY.fa"),
  genes = 
    c("/Users/sarahcarey/Desktop/hop_genomes/annotations/Humulus_lupulus_lupulus_21110M_HAP2_v1.2.primary.X.gff3",
      "/Users/sarahcarey/Desktop/hop_genomes/annotations/Humulus_lupulus_lupulus_21110M_HAP1_v1.3.primary.Y.gff3"),
  reps = 
    c("/Users/sarahcarey/Desktop/hop_genomes/annotations/repeats/EDTA/Humulus_lupulus_lupulus_21110M_HAP2_v1.0.fa.mod.EDTA.TEanno.X.gff3",
      "/Users/sarahcarey/Desktop/hop_genomes/annotations/repeats/EDTA/Humulus_lupulus_lupulus_21110M_HAP1_v1.2.fa.mod.EDTA.TEanno.Y.gff3"))

Ymers <- Biostrings::readDNAStringSet("/Users/sarahcarey/Desktop/hop_genomes/final_assemblies/lupulus_maleOnly_7.21mers.list.txt.fasta")
Ymers <- as.character(Ymers)

test <- plot_2genomes(
  wd = "/Users/sarahcarey/Desktop/hop_genomes/genespace",
  genomeIDs <- c("21110M_X", "21110M_Y"),
  faFiles = filepaths$assemblies,
  geneGffFiles = filepaths$genes,
  repeatGffFiles = filepaths$reps,
  nCores = 12, repeatClassColumnName = "type",
  repeatGrep1 = "Gypsy",
  repeatGrep2 = "Copia", returnSourceData = TRUE,
  overwrite = TRUE,
  forceCleanWindows = FALSE, kmers = Ymers)




##### Figure 2a, synteny between hops, cannabis, mulberries, and figs #####

library("GENESPACE")
#v1.3.1

wd <- "/Users/sarahcarey/Desktop/hop_genomes/Cannabaceae_genespace_June2025"
genomeRepo <-"/Users/sarahcarey/Desktop/genespace/genomes"
path2mcscanx <- "/Users/sarahcarey/Desktop/packages/MCScanX-master"

parse_annotations(rawGenomeRepo=genomeRepo,
                  genomeDirs=c("Cannabis_OttoII_HAP2"),
                  genomeIDs = c("Cannabis_OttoII_HAP2"),
                  gffString = "gff3",
                  faString = "fa",
                  genespaceWd=wd,
                  troubleShoot = FALSE,
                  headerEntryIndex = 1,
                  overwrite = F,
                  headerSep=" ",
                  gffIdColumn = "ID")

parse_annotations(rawGenomeRepo=genomeRepo,
                  genomeDirs=c("lupulus_hap2"),
                  genomeIDs = c("Humulus_lupulus_HAP2"),
                  gffString = "gff3",
                  faString = "fa",
                  genespaceWd=wd,
                  troubleShoot = FALSE,
                  headerEntryIndex = 1,
                  overwrite = F,
                  headerSep=" ",
                  gffIdColumn = "ID")


parse_annotations(rawGenomeRepo=genomeRepo,
                  genomeDirs="MorusAlba",
                  genomeIDs = "Morus_alba",
                  gffString = "gff",
                  faString = "faa",
                  genespaceWd=wd,
                  troubleShoot = FALSE,
                  headerEntryIndex = 1,
                  overwrite = F,
                  headerSep="\t",
                  gffIdColumn = "Protein_Accession")

parse_annotations(rawGenomeRepo=genomeRepo,
                  genomeDirs="FicusHispida",
                  genomeIDs = "Ficus_hispida",
                  gffString = "gff",
                  faString = "faa",
                  genespaceWd=wd,
                  troubleShoot = FALSE,
                  headerEntryIndex = 1,
                  overwrite = F,
                  headerSep="\t",
                  gffIdColumn = "Protein_Accession")


gpar <- init_genespace(genomeIDs = c("Humulus_lupulus_HAP2","Cannabis_OttoII_HAP2","Morus_alba", "Ficus_hispida"), 
                       ploidy = c(1,1,1,1),
                       wd = wd, nCores = 12, path2mcscanx = path2mcscanx)


roi <- data.frame(genome = c("Humulus_lupulus_HAP2","Ficus_hispida", "Morus_alba"), 
                  chr = c("ChrX","GWHALOG00000001", "GWHDOOQ00000001"), 
                  color = c("seagreen","salmon","#FF0000"), start =c(1,1,1), 
                  end = c(300000000,2000000000,20000000000))


ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))

pdf("/Users/sarahcarey/Desktop/hop_genomes/Cannabaceae_genespace_June2025/swPlot_combRaw2.pdf",
    height = 6, width = 8)
ripDat <- plot_riparian(gpar, refGenome = "Humulus_lupulus_HAP2",
                        chrLabFontSize = 7,
                        genomeIDs = c("Morus_alba",
                                      "Ficus_hispida",
                                      "Cannabis_OttoII_HAP2",
                                      "Humulus_lupulus_HAP2"), 
                        highlightBed = roi,
                        addThemes = ggthemes, 
                        labelTheseGenomes = c("Humulus_lupulus_HAP2","Cannabis_OttoII_HAP2"),
                        chrFill = "black", minChrLen2plot = 200)
dev.off()


##### Figure 2b, synteny between hops, Ks of gametologs #####

### GENESPACE between hop and hemp

library("GENESPACE")
#v1.3.1

filepaths <- list(
  assemblies = c(
    "/Users/sarahcarey/Desktop/hop_genomes/final_assem/lupulus_ChrX.rc.fasta",
    "/Users/sarahcarey/Desktop/hemp_genomes/final_assem/otto_ChrX.fasta"),
  genes = c(
    "/Users/sarahcarey/Desktop/hop_genomes/annotations/Humulus_lupulus_lupulus_21110M_HAP2_v1.2.X.gff3",
    "/Users/sarahcarey/Desktop/hemp_genomes/OIIb.primary_high_confidence.ChrRenamed.X.gff3"),
  reps = c(
    "/Users/sarahcarey/Desktop/hop_genomes/Humulus_lupulus_lupulus_21110M_HAP2_v1.0.fa.mod.EDTA.TEanno.X.gff3",
    "/Users/sarahcarey/Desktop/hemp_genomes/Cannabis_sativa_OttoII_HAP2_v1.0.fa.mod.EDTA.TEanno.X.gff3"))


test <- plot_2genomes(
  wd = "/Users/sarahcarey/Desktop/hop_genomes/final_assem/genespace",
  genomeIDs <- c("lupulus_X", "otto_X"),
  faFiles = filepaths$assemblies,
  geneGffFiles = filepaths$genes, 
  repeatGffFiles = filepaths$reps,
  nCores = 12, repeatClassColumnName = "type",
  repeatGrep1 = "Gypsy",
  repeatGrep2 = "Copia", returnSourceData = TRUE,  
  overwrite = TRUE, minChrSize = 5000000, windowSize=1000, 
  forceCleanWindows = TRUE, mm2mode="default", asmPreset="asm10")



## Ks plots

library(karyoploteR)
#v1.32.0

## read in data
sex_chroms <- read.table("figure_data/Cannabaceae.sexChroms.lengths.karyoploter.txt", header=TRUE)
sex_chroms.GR <- toGRanges(data.frame(chr=sex_chroms$Chromosome, start=sex_chroms$ChromStart, 
                                      end=sex_chroms$ChromEnd))

coverage_data_Ks <- read.csv("figure_data/lupulus_otto_Ks.csv", header=TRUE)
coverage_data_Ks.GR <- toGRanges(data.frame(chr=coverage_data_Ks$Chromosome, 
                                                    start=coverage_data_Ks$Start, end=coverage_data_Ks$Stop))

pp <- getDefaultPlotParams(plot.type = 1)
#pp$data2height <- 0
pp$data1outmargin <- 100
#pp$data2outmargin <- 0


pdf("figure_output/Ks_plot.pdf")

kp <- plotKaryotype(sex_chroms.GR, pin=8, plot.type = 1, labels.plotter = NULL, plot.params = pp, 
                    chromosomes = c("lupulus_ChrX", "otto_ChrX"))
kpAddChromosomeNames(kp, xoffset=-0.01, cex=0.75)
kpAddBaseNumbers(kp, tick.dist = 25000000, tick.len = 10, tick.col="black", cex=0.75,
                 minor.tick.dist = 10000000, minor.tick.len = 10, minor.tick.col = "black", add.units = F)

kp <-kpPoints(kp, data=coverage_data_Ks.GR, data.panel = 1, 
               y=coverage_data_Ks$Ks, col=transparent("black", amount=0.5), ymin=0,
              ymax=0.8)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, cex=0.75, numticks=3,
       labels=c("0", "0.4", "0.8"), side=1)

dev.off()


##### Figure 2c, synteny between hops, Ks of gametologs #####
## analysis by Phil Bentz



##### Figure 3a, repeat Kimura substitution plots #####


library("ggplot2")
#v3.5.2

library("reshape2")
#v1.4.4

library("viridis")
#v0.6.5

library("ggpubr")
#v0.6.1

KimuraDistance_LuY <- read.csv("figure_data/Humulus_lupulus_lupulus_21110M_HAP1_v1.2.ChrY.tbl",sep=" ")

KimuraDistance_LuY$OtherRepeat <- (KimuraDistance_LuY$DNA.DTA+ 
                                     KimuraDistance_LuY$DNA.DTC+ KimuraDistance_LuY$DNA.DTH+ 
                                     KimuraDistance_LuY$DNA.DTM+ KimuraDistance_LuY$DNA.DTT+ 
                                     KimuraDistance_LuY$DNA.Helitron+ KimuraDistance_LuY$LINE.unknown+ 
                                     KimuraDistance_LuY$LTR.unknown+ KimuraDistance_LuY$MITE.DTA+ 
                                     KimuraDistance_LuY$MITE.DTC+ KimuraDistance_LuY$MITE.DTH+ 
                                     KimuraDistance_LuY$MITE.DTM+ KimuraDistance_LuY$MITE.DTT+ 
                                     KimuraDistance_LuY$Simple_repeat+ KimuraDistance_LuY$Unknown)
KimuraDistance_LuY_combined <- as.data.frame(cbind(KimuraDistance_LuY$Div, 
                                     KimuraDistance_LuY$LTR.Gypsy, KimuraDistance_LuY$LTR.Copia,
                                     KimuraDistance_LuY$OtherRepeat))
colnames(KimuraDistance_LuY_combined) <- c("Div", "Ty-3", "Ty-1", "Other repeat")

KimuraDistance_LuX <- read.csv("figure_data/Humulus_lupulus_lupulus_21110M_HAP2_v1.0.ChrX.tbl",sep=" ")

KimuraDistance_LuX$OtherRepeat <- (KimuraDistance_LuX$DNA.DTA+ 
                                     KimuraDistance_LuX$DNA.DTC+ KimuraDistance_LuX$DNA.DTH+ 
                                     KimuraDistance_LuX$DNA.DTM+ KimuraDistance_LuX$DNA.DTT+ 
                                     KimuraDistance_LuX$DNA.Helitron+ KimuraDistance_LuX$LINE.unknown+ 
                                     KimuraDistance_LuX$LTR.unknown+ KimuraDistance_LuX$MITE.DTA+ 
                                     KimuraDistance_LuX$MITE.DTC+ KimuraDistance_LuX$MITE.DTH+ 
                                     KimuraDistance_LuX$MITE.DTM+ KimuraDistance_LuX$MITE.DTT+ 
                                     KimuraDistance_LuX$Simple_repeat+ KimuraDistance_LuX$Unknown)
KimuraDistance_LuX_combined <- as.data.frame(cbind(KimuraDistance_LuX$Div, 
                                                   KimuraDistance_LuX$LTR.Gypsy, KimuraDistance_LuX$LTR.Copia,
                                                   KimuraDistance_LuX$OtherRepeat))
colnames(KimuraDistance_LuX_combined) <- c("Div", "Ty-3", "Ty-1", "Other repeat")

KimuraDistance_NeY <- read.csv("figure_data/Humulus_lupulus_neomexicanus_MN-586_HAP1_v1.0.ChrY.tbl",sep=" ")

KimuraDistance_NeY$OtherRepeat <- (KimuraDistance_NeY$DNA.DTA+ 
                                     KimuraDistance_NeY$DNA.DTC+ KimuraDistance_NeY$DNA.DTH+ 
                                     KimuraDistance_NeY$DNA.DTM+ KimuraDistance_NeY$DNA.DTT+ 
                                     KimuraDistance_NeY$DNA.Helitron+ KimuraDistance_NeY$LINE.unknown+ 
                                     KimuraDistance_NeY$LTR.unknown+ KimuraDistance_NeY$MITE.DTA+ 
                                     KimuraDistance_NeY$MITE.DTC+ KimuraDistance_NeY$MITE.DTH+ 
                                     KimuraDistance_NeY$MITE.DTM+ KimuraDistance_NeY$MITE.DTT+ 
                                     KimuraDistance_NeY$Simple_repeat+ KimuraDistance_NeY$Unknown)
KimuraDistance_NeY_combined <- as.data.frame(cbind(KimuraDistance_NeY$Div, 
                                                   KimuraDistance_NeY$LTR.Gypsy, KimuraDistance_NeY$LTR.Copia,
                                                   KimuraDistance_NeY$OtherRepeat))
colnames(KimuraDistance_NeY_combined) <- c("Div", "Ty-3", "Ty-1", "Other repeat")

KimuraDistance_NeX <- read.csv("figure_data/Humulus_lupulus_neomexicanus_MN-586_HAP2_v1.0.ChrX.tbl",sep=" ")

KimuraDistance_NeX$OtherRepeat <- (KimuraDistance_NeX$DNA.DTA+ 
                                     KimuraDistance_NeX$DNA.DTC+ KimuraDistance_NeX$DNA.DTH+ 
                                     KimuraDistance_NeX$DNA.DTM+ KimuraDistance_NeX$DNA.DTT+ 
                                     KimuraDistance_NeX$DNA.Helitron+ KimuraDistance_NeX$LINE.unknown+ 
                                     KimuraDistance_NeX$LTR.unknown+ KimuraDistance_NeX$MITE.DTA+ 
                                     KimuraDistance_NeX$MITE.DTC+ KimuraDistance_NeX$MITE.DTH+ 
                                     KimuraDistance_NeX$MITE.DTM+ KimuraDistance_NeX$MITE.DTT+ 
                                     KimuraDistance_NeX$Simple_repeat+ KimuraDistance_NeX$Unknown)
KimuraDistance_NeX_combined <- as.data.frame(cbind(KimuraDistance_NeX$Div, 
                                                   KimuraDistance_NeX$LTR.Gypsy, KimuraDistance_NeX$LTR.Copia,
                                                   KimuraDistance_NeX$OtherRepeat))
colnames(KimuraDistance_NeX_combined) <- c("Div", "Ty-3", "Ty-1", "Other repeat")




#add the genome size in bp
#lupulus Y
genomes_size_LuY=168676604
#lupulus X
genomes_size_LuX=268239420

#neomexicanus Y
genomes_size_NeY=246096435
#neomexicanus X
genomes_size_NeX=267359615

kd_melt_LuY = melt(KimuraDistance_LuY_combined,id="Div")
kd_melt_LuY$norm = kd_melt_LuY$value/genomes_size_LuY * 100
kd_melt_NeY = melt(KimuraDistance_NeY_combined,id="Div")
kd_melt_NeY$norm = kd_melt_NeY$value/genomes_size_NeY * 100


kd_melt_LuX = melt(KimuraDistance_LuX_combined,id="Div")
kd_melt_LuX$norm = kd_melt_LuX$value/genomes_size_LuX * 100
kd_melt_NeX = melt(KimuraDistance_NeX_combined,id="Div")
kd_melt_NeX$norm = kd_melt_NeX$value/genomes_size_NeX * 100


# using Time = Kimura Div / (2 * mu)
kd_melt_LuY$time = ((kd_melt_LuY$Div/100)/(2 * 1e-9)/1e6)
kd_melt_LuX$time = ((kd_melt_LuX$Div/100)/(2 * 1e-9)/1e6)
kd_melt_NeY$time = ((kd_melt_NeY$Div/100)/(2 * 1e-9)/1e6)
kd_melt_NeX$time = ((kd_melt_NeX$Div/100)/(2 * 1e-9)/1e6)



### make the plots

dup_labels <- function(x) {kd_melt_LuY$time[match(x, kd_melt_LuY$Div)]}

LuY <- ggplot(kd_melt_LuY, aes(fill=variable, y=norm, x=Div)) + 
  geom_bar(position="stack", stat="identity", width = 1) +
  theme_classic() +
  xlab("") +
  ylab("Percent of the genome") + 
  labs(fill = "") +
  coord_cartesian(xlim = c(0, 50)) +
  theme(axis.text=element_text(size=8),axis.title =element_text(size=10)) +
  ylim(0,6) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3), 
                     sec.axis = dup_axis(labels = dup_labels, name = "")) +
  geom_vline(xintercept=9.6) +
  scale_fill_manual(values = c("#052d61", "#1e90ff", "#87ceeb")) + 
  theme(legend.position="none")


NeY <- ggplot(kd_melt_NeY, aes(fill=variable, y=norm, x=Div))  + 
  geom_bar(position="stack", stat="identity", width = 1) +
  theme_classic() +
  xlab("") +
  ylab("") + 
  labs(fill = "") +
  coord_cartesian(xlim = c(0, 50)) +
  theme(axis.text=element_text(size=8),axis.title =element_text(size=10)) +
  ylim(0,6) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3), 
                     sec.axis = dup_axis(labels = dup_labels, name = "")) +
  geom_vline(xintercept=9.6) +
  scale_fill_manual(values = c("#052d61", "#1e90ff", "#87ceeb")) + 
  theme(legend.position="none")


## Xs
LuX <- ggplot(kd_melt_LuX, aes(fill=variable, y=norm, x=Div))  + 
  geom_bar(position="stack", stat="identity", width = 1) +
  theme_classic() +
  xlab("") +
  ylab("") + 
  labs(fill = "") +
  coord_cartesian(xlim = c(0, 50)) +
  theme(axis.text=element_text(size=8),axis.title =element_text(size=10)) +
  ylim(0,6) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3), 
                     sec.axis = dup_axis(labels = dup_labels, name = "")) +
  geom_vline(xintercept=9.6) +
  scale_fill_manual(values = c("#052d61", "#1e90ff", "#87ceeb")) + 
  theme(legend.position="none")

NeX <- ggplot(kd_melt_NeX, aes(fill=variable, y=norm, x=Div))  + 
  geom_bar(position="stack", stat="identity", width = 1) +
  theme_classic() +
  xlab("") +
  ylab("") + 
  labs(fill = "") +
  coord_cartesian(xlim = c(0, 50)) +
  theme(axis.text=element_text(size=8),axis.title =element_text(size=10)) +
  ylim(0,6) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3), 
                     sec.axis = dup_axis(labels = dup_labels, name = "")) +
  geom_vline(xintercept=9.6) +
  scale_fill_manual(values = c("#052d61", "#1e90ff", "#87ceeb")) + 
  theme(legend.position="none")

### arrange the plots

kimura_combined <- ggarrange(LuY, NeY, 
                             LuX, NeX, 
                             labels=c(""),
                             ncol=2, nrow=2)
kimura_combined

ggsave(kimura_combined, filename = "figure_output/kimura_combined_reduced.pdf", device = pdf, dpi = 300,
       width = 3.5, height = 3, units = "in")



##### Figure 3b, dN/dS boxplots #####

library("ggplot2")
#v3.5.2

library("ggsignif")
#v0.6.4


## load in data
paml_result <- read.csv("figure_data/paml.csv", header= TRUE)

## convert to numeric
paml_result$dN.dS <- as.numeric(paml_result$dN.dS)

## wilcox test, to assess significance (added to plot using geom_signif)
paml_result_cannabis <- subset(paml_result, paml_result$species=="cannabis")
paml_result_hops <- subset(paml_result, paml_result$species=="hop")

pairwise.wilcox.test(paml_result_cannabis$dN.dS,paml_result_cannabis$bin_adv,
                     p.adjust.method = "BH")
#Autosomal      HXR    
#HXR 8.3e-08     -      
#PAR 0.05       8.2e-08


pairwise.wilcox.test(paml_result_hops$dN.dS,paml_result_hops$bin_adv,
                     p.adjust.method = "BH")
#Autosomal    HXR    
#HXR 0.00018   -      
#PAR 0.98615   0.02055

## plot
p <- ggplot(paml_result, aes(x=species, y=dN.dS)) +
  ylab("dN/dS") +
  theme(axis.title.x = element_text(size=10)) +
  theme(axis.title.y = element_text(size=10)) +
  theme(axis.text.x = element_text(size=5),
        axis.text.y = element_text(size=5)) +
  theme(panel.background = element_rect(fill="white")) +
  guides(fill=FALSE) + 
  theme(axis.text.x= element_text(size=5)) +
  geom_point(position=position_jitterdodge(jitter.width = 0.5, jitter.height = 0), 
             aes(group=bin_adv, color=factor(bin_adv)), alpha=0.25, shape=19, size=1) + 
  theme(legend.position="right") +
  scale_colour_manual("Genomic location", values=c("#FFC520", "#008656", "#58B947")) + 
  geom_boxplot(aes(colour = bin_adv), fill="white", alpha=0.25, lwd=1) + 
  geom_signif(
    y_position = c(2.5, 2.75, 3, 2.5, 2.75, 2.9), xmin = c(0.75, 1, 0.75, 1.75, 2, 1.75), 
    xmax = c(1, 1.25, 1.25, 2, 2.25, 2.25),
    annotation = c("***", "***", "*", "***", "*", "ns"), tip_length = 0.005, textsize = 3, size = 1)
p

ggsave(p, filename = "figure_output/fasterX.pdf", device = pdf, dpi = 300,
       width = 4, height = 3.5, units = "in")



##### Figure 3c, Fst outliers, chemotypes #####

library(karyoploteR)
#v1.32.0

## read in data
otto_chroms <- read.table("figure_data/OttoII_HAP2_karyoploter.txt", header=TRUE)
otto_chroms.GR <- toGRanges(data.frame(chr=otto_chroms$Chromosome, start=otto_chroms$ChromStart, 
                                      end=otto_chroms$ChromEnd))

coverage_data_Fst <- read.csv("figure_data/cannabis_FstOutliers.csv", header=TRUE)
coverage_data_Fst.GR <- toGRanges(data.frame(chr=coverage_data_Fst$chromosome, 
                                            start=coverage_data_Fst$window_pos_1, end=coverage_data_Fst$window_pos_2))

pp <- getDefaultPlotParams(plot.type = 4)
pp$ideogramlateralmargin <- 0.01
pp$leftmargin <- 0.1


pdf("figure_output/Fst_chemotypes_plot_ChrX.pdf", width = 7.5, height = 2.5)

kp <- plotKaryotype(otto_chroms.GR, pin=8, plot.type = 4, labels.plotter = NULL, plot.params = pp, 
                    chromosomes = "ChrX")

kpAddChromosomeNames(kp, yoffset=0)
kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 5, tick.col="black", cex=0.5,
                 minor.tick.dist = 5000000, minor.tick.len = 5, minor.tick.col = "black", add.units = F)

kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, cex=0.75, numticks=3,
       labels=c("0", "0.25", "0.5"), side=1)

markers_X <- data.frame(chr=rep("ChrX",9), 
                        pos=(c(29341994,
                               29349058,
                               29388947,
                               29463831,
                               29613142,
                               20539916,
                               26228604,
                               24464515,
                               19996372)), 
                        labels=paste0(c("PT",
                                        "PT",
                                        "PT",
                                        "PT",
                                        "PT",
                                        "PT",
                                        "FT",
                                        "GI",
                                        "CO")))

kpPlotMarkers(kp, chr=markers_X$chr, x=markers_X$pos, labels=markers_X$labels, 
              cex=0.75, text.orientation = "horizontal", r1=1, clipping = FALSE)


kp <-kpPoints(kp, data=coverage_data_Fst.GR, data.panel = 1, 
              y=coverage_data_Fst$fst_chemotypes_femalesOnly, 
              col=colByCategory(coverage_data_Fst$fst_chemotypes_femalesOnly_outlier, c("gray75","black")), 
              ymin=0,ymax=0.5)


kp <-kpLines(kp, data=coverage_data_Fst.GR, data.panel = 1, 
              y=.206, 
              col= "black", 
              ymin=0,ymax=0.5, lwd=2)

dev.off()




##### Figure 4, female vs monoecious analyses #####


library(karyoploteR)
#v1.32.0

## read in data
otto_chroms <- read.table("figure_data/OttoII_HAP2_karyoploter.txt", header=TRUE)
otto_chroms.GR <- toGRanges(data.frame(chr=otto_chroms$Chromosome, start=otto_chroms$ChromStart, 
                                       end=otto_chroms$ChromEnd))

coverage_data_Fst <- read.csv("figure_data/cannabis_FstOutliers.csv", header=TRUE)
coverage_data_Fst.GR <- toGRanges(data.frame(chr=coverage_data_Fst$chromosome, 
                                             start=coverage_data_Fst$window_pos_1, end=coverage_data_Fst$window_pos_2))

coverage_data_hits <- read.csv("figure_data/Khufu_results.csv", header=TRUE)
coverage_hits.GR <- toGRanges(data.frame(chr=coverage_data_hits$Chromosome, start=coverage_data_hits$Position, end=coverage_data_hits$Position))
coverage_data_hits_subset_0to4 <- subset(coverage_data_hits, coverage_data_hits$Contrast=="0to4")
coverage_data_hits_subset_0to3 <- subset(coverage_data_hits, coverage_data_hits$Contrast=="0to3")
coverage_data_hits_subset_0to123 <- subset(coverage_data_hits, coverage_data_hits$Contrast=="0toAllMonoecious")


gene_expr <- read.csv("figure_data/geneExpr_meristemSexes.csv", header=TRUE)
gene_expr.GR <- toGRanges(data.frame(chr=gene_expr$Chromosome, 
                                     start=gene_expr$Start, end=gene_expr$Stop))

filtered.dm.genes <- na.omit(gene_expr)
log.pval <- as.data.frame(-log10(filtered.dm.genes$padj))
filtered.dm.genes <- cbind(filtered.dm.genes, log.pval)
filtered.dm.genes

sign.genes <- filtered.dm.genes[filtered.dm.genes$padj < 0.05,]


pp <- getDefaultPlotParams(plot.type = 4)
pp$ideogramlateralmargin <- 0.01
pp$leftmargin <- 0.1


#### panel a

pdf("figure_output/Fst_Khufu_FstMonoecy_AllChroms_plot.pdf", width = 7.5, height = 3.5)

kp <- plotKaryotype(otto_chroms.GR, pin=8, plot.type = 4, labels.plotter = NULL, plot.params = pp)

kpAddChromosomeNames(kp, yoffset=-5)
kpAddBaseNumbers(kp, tick.dist = 25000000, tick.len = 5, tick.col="black", cex=0.5,
                 minor.tick.dist = 5000000, minor.tick.len = 5, minor.tick.col = "black", add.units = F)

kp <-kpPoints(kp, data=coverage_hits.GR, data.panel = 1, col="black",
              y=0.95, x=coverage_data_hits_subset_0to123$Position, ymin=0, ymax=1, r0=0.8, r1=0.9, lwd=2, cex=1)


kp <-kpPoints(kp, data=coverage_data_Fst.GR, data.panel = 1, 
              y=coverage_data_Fst$fst_FvsMono, 
              col=colByCategory(coverage_data_Fst$fst_FvsMono_outlier, c("gray50","black")), 
              ymin=0, ymax=0.8, r0=0, r1=0.65)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, cex=0.75, numticks=3,
       labels=c("0", "0.4", "0.8"), side=1, r0=0, r1=0.65)

kp <-kpLines(kp, data=coverage_data_Fst.GR, data.panel = 1, 
             y=.18, 
             col= "black", 
             ymin=0,ymax=0.5, lwd=2)

dev.off()



#### panel b

pdf("figure_output/Fst_Khufu_FstMonoecy_ChrX_plot.pdf", width = 7.5, height = 3.5)

zoom.region <- toGRanges(data.frame("ChrX", 76000000, 86000000))
kp <- plotKaryotype(otto_chroms.GR, pin=8, plot.type = 4, labels.plotter = NULL, plot.params = pp,
                    chromosomes = "ChrX", zoom = zoom.region)

kpAddChromosomeNames(kp, yoffset=-1)
kpAddBaseNumbers(kp, tick.dist = 1000000, tick.len = 5, tick.col="black", cex=0.5,
                 minor.tick.dist = 5000000, minor.tick.len = 5, minor.tick.col = "black", add.units = F)

kpRect(kp, chr="ChrX", x0=77637719, x1=81741939, y0=0, y1=1, col="gray90", data.panel=1, border=NA)

#kpRect(kp, chr="ChrX", x0=79228107, x1=81743844, y0=0, y1=1, col="gray50", data.panel=1, border=NA)


markers_X <- data.frame(chr=rep("ChrX", 9), 
                        pos=(c(79491408,
                               69758278,
                               85646292,
                               79420542,
                               49744028,
                               71576472,
                               74717410,
                               80781598,
                               80760243)), 
                        labels=paste0(c("ACS",
                                        "FD",
                                        "ETR",
                                        "ACO",
                                        "EIL",
                                        "EBF",
                                        "CTR",
                                        "REM",
                                        "KAN")))

kpPlotMarkers(kp, chr=markers_X$chr, x=markers_X$pos, labels=markers_X$labels, 
              cex=0.75, text.orientation = "horizontal", r1=1, clipping = FALSE)

kp <-kpPoints(kp, data=coverage_hits.GR, data.panel = 1, col="black",
              y=0.95, x=coverage_data_hits_subset_0to123$Position, ymin=0, ymax=1, r0=0.8, r1=0.9, lwd=2, cex=1)


kp <-kpPoints(kp, data=coverage_data_Fst.GR, data.panel = 1, 
              y=coverage_data_Fst$fst_FvsMono, 
              col=colByCategory(coverage_data_Fst$fst_FvsMono_outlier, c("gray50","black")), 
              ymin=0, ymax=0.8, r0=0.35, r1=0.65)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, cex=0.75, numticks=3,
       labels=c("0", "0.4", "0.8"), side=1, r0=0.35, r1=0.65)

kp <-kpLines(kp, data=coverage_data_Fst.GR, data.panel = 1, 
             y=.26, 
             col= "black", 
             ymin=0,ymax=0.5, lwd=2)

#fc.ymax <- ceiling(max(abs(range(sign.genes$log2FoldChange))))
#fc.ymin <- -fc.ymax

col.over <- "#008656"
col.under <-"#FFC520"


sign.col <- rep(col.over, 3642)
sign.col[sign.genes$log2FoldChange<0] <- col.under
#cex.val <- sqrt(sign.genes$`-log10(filtered.dm.genes$padj)`)/2
#points.top <- 0.8

kpPoints(kp, chr= sign.genes$Chromosome, x = sign.genes$Start, 
         y=sign.genes$log2FoldChange, ymax=15, ymin=-15, 
         cex=1, col=sign.col, r1=0.3, r0=0)
kpAxis(kp, ymax=kp$latest.plot$computed.values$max.density, cex=0.75, numticks=3,
       labels=c("-15", "0", "15"), side=1,  r1=0.3, r0=0)


dev.off()

################# SUPPLEMENTAL FIGURES #################



#### Change points ####

library("mcp") 
#v0.3.4
library("rjags")
#v4-17
library("PMCMRplus")
#v1.9.12
library("ggplot2")
#v3.5.2

Ks <- read.csv("figure_data/lupulus_otto_Ks.csv", header=TRUE)
Ks_hemp <- subset(Ks, Chromosome == "otto_ChrX")
#Ks_hemp_SDR <- subset(Ks_hemp, Compartment == "SDR")
Ks_hop <- subset(Ks, Chromosome == "lupulus_ChrX")
#Ks_hop_SDR <- subset(Ks_hop, Compartment == "SDR")


## changepoint using Otto II
Ks_hemp_df <- as.data.frame(Ks_hemp$Start)
Ks_hemp_df <- as.data.frame(cbind(Ks_hemp_df,Ks_hemp$Ks))
colnames(Ks_hemp_df) <- c("x", "y")

## checking for outliers
out_hemp <- gesdTest(Ks_hemp_df$y, maxr=20)
summary(out_hemp)
boxplot(Ks_hemp_df$y)

## removing significant outliers
Ks_hemp_df <- Ks_hemp_df[-119,]

#one change (two plateaus)
model = list(y ~ 1, ~ 1)

#two changes (three plateaus)
model = list(y ~ 1, ~ 1, ~ 1)

#three changes (four plateaus)
model = list(y ~ 1, ~ 1, ~ 1, ~ 1)

fit_hemp = mcp(model, Ks_hemp_df, par_x = "x", sample = "both", iter = 10000, adapt=10000, chains = 3)

png("figure_output/cannabis_changepoint_Ks_3.png", width = 8, height = 4, units = 'in', res = 500)
plot(fit_hemp , q_fit = TRUE, lines=50) + ggtitle("Posterior fit") + 
  xlab("Position (Mb)") + ylab("Ks") +
  scale_x_continuous(labels = function(x)x/1000000, n.breaks = 20)
dev.off()

out_hemp <- summary(fit_hemp)
write.csv(out_hemp, "cannabis_changepoint_Ks.csv")

## changepoint using hops
Ks_hop_df <- as.data.frame(Ks_hop$Start)
Ks_hop_df <- as.data.frame(cbind(Ks_hop_df,Ks_hop$Ks))
colnames(Ks_hop_df) <- c("x", "y")

## checking for outliers
out_hop <- gesdTest(Ks_hop_df$y, maxr=5)
summary(out_hop)
boxplot(Ks_hop_df$y)

#model = list(y ~ 1)
model = list(y ~ 1, ~ 1)
model = list(y ~ 1, ~ 1, ~ 1, ~ 1)

fit_hop = mcp(model, Ks_hop_df, par_x = "x", sample = "both", iter = 10000, adapt=10000, chains = 3)

png("figure_output/hop_changepoint_Ks_3.png", width = 8, height = 4, units = 'in', res = 500)
plot(fit_hop , q_fit = TRUE, lines=50) + ggtitle("Posterior fit") + 
  xlab("Position (Mb)") + ylab("Ks") +
  scale_x_continuous(labels = function(x)x/1000000, n.breaks = 20)
dev.off()


out_hop <- summary(fit_hop)
write.csv(out_hop, "cannabis_changepoint_Ks.csv")


##### gene expression heatmap #####

library("DESeq2")
#v1.46.0
library("pheatmap")
#v1.10.12
library("ggplot2")
#v3.5.2


countData <- as.matrix(read.csv("figure_data/heatmap_gene_count_matrix.csv", row.names = "gene_id"))
colData <- read.csv("figure_data/heatmap_pheno.csv", sep=",", header=TRUE) 


dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~Sex+Tissue)
dds <- DESeq(dds)

ntd <- normTransform(dds)

## ACS/ETR/ACOs/REM/KAN
select <- c("OIIb.chrX.v1.g435690",
            "OIIb.chr4.v1.g213660",
            "OIIb.chr5.v1.g223160",
            "OIIb.chr1.v1.g077020",  
            "OIIb.chr1.v1.g091300",
            "OIIb.chr5.v1.g252740",
            "OIIb.chr3.v1.g172660",
            "OIIb.chr5.v1.g246710",
            "OIIb.chrX.v1.g440030",
            "OIIb.chr7.v1.g309190",
            "OIIb.chrX.v1.g435640",
            "OIIb.chr1.v1.g077290",
            "OIIb.chr1.v1.g100060", 
            "OIIb.chr8.v1.g343400",
            "OIIb.chrX.v1.g436840",
            "OIIb.chrX.v1.g436820",
            "OIIb.chrX.v1.g436780")


df <- as.data.frame(colData(dds)[,c("Sex","Tissue")])
rownames(df) <- colnames(dds)
colnames(df) <- c("Sex","Tissue")

pheatmap(assay(ntd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, cex=0.75)

df <- as.data.frame(colData(dds)[,c("Sex","Tissue")])
rownames(df) <- colnames(dds)
colnames(df) <- c("Sex","Tissue")

heatmap <- pheatmap(assay(ntd)[select,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df, cex=0.75)


ggsave(heatmap, filename = "figure_output/heatmap.pdf", device = pdf, dpi = 300,
       width = 12, height = 6, units = "in")


##### Khufu Y chromosome test #####

## panel a by Zack Myers

## panel b

library("ggplot2")
#v3.5.2

isolate_sex <- read.csv("Figure_data/hemp_NWG_Khufu_sdrCov.csv", header=TRUE)

isolate_sex$Isolate <- factor(isolate_sex$Genotype, 
                              levels = isolate_sex$Genotype[order(isolate_sex$ChrY_cov, 
                                                                  decreasing = TRUE)])
p <- ggplot(isolate_sex, aes(x=Genotype, y=ChrY_cov)) + 
  geom_point(alpha=0.75, size=3,(aes(col=Sex, shape=Sex))) +
  xlab("Isolate") +
  ylab("Coverage") +
  theme(axis.title.x = element_text(size=10)) +
  theme(axis.title.y = element_text(size=10))+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=5),
        axis.text.y = element_text(size=15)) +
  theme(plot.margin=unit(c(1,1,1,1),"cm")) + 
  theme(panel.background = element_rect(color="black")) + 
  theme(axis.text.x=element_blank())
p

ggsave("figure_output/hemp_isolate_SDR_cov.png", p, units="in", width=8, height=4, dpi=300,
       device="png")




######## Figs. S3-4, contig maps #########

library("GENESPACE")
#v1.3.1
library("Biostrings")
#v2.74

#for hop

faf <- c("/Users/sarahcarey/Desktop/hop_genomes/final_assem/Humulus_lupulus_lupulus_21110M_HAP1_v1.2.fa",
         "/Users/sarahcarey/Desktop/hop_genomes/final_assem/Humulus_lupulus_lupuloides_MN-1421_HAP1_v1.2.fa",
         "/Users/sarahcarey/Desktop/hop_genomes/final_assem/Humulus_lupulus_neomexicanus_MN-586_HAP1_v1.0.fa",
         "/Users/sarahcarey/Desktop/hop_genomes/final_assem/Humulus_lupulus_lupulus_21110M_HAP2_v1.0.fa",
         "/Users/sarahcarey/Desktop/hop_genomes/final_assem/Humulus_lupulus_lupuloides_MN-1421_HAP2_v1.0.fa",
         "/Users/sarahcarey/Desktop/hop_genomes/final_assem/Humulus_lupulus_neomexicanus_MN-586_HAP2_v1.0.fa")

names(faf) <- c("HlLuH1", "HlLoH1", "HlNeH1", "HlLuH2", "HlLoH2", "HlNeH2") 

ssl <- sapply(faf, USE.NAMES = T, simplify = F, function(x){
  y <- readDNAStringSet(x)
  return(y)
})

cgtl <- sapply(names(ssl), simplify = F, USE.NAMES = T, function(i){
  print(i)
  find_contigsGapsTelos(
    ssl[[i]],
    teloKmers = c("TTTAGGG"),
    minChrSize = 1e7,
    maxDistBtwTelo = 100,
    minTeloSize = 400,
    minTeloDens = 0.75,
    maxDist2end = 50000)
})
p1 <- plot_contigs(
  cgtList = cgtl,
  palette = viridis::turbo, nColors = 5)

pdf("Hop_contigs_telos.pdf", width = 10, height = 8)
p1 <- plot_contigs(
  cgtList = cgtl,
  palette = viridis::turbo, nColors = 5)
dev.off()

#for cannabis

faf <- c("/Users/sarahcarey/Desktop/hemp_genomes/final_assem/Cannabis_sativa_Carmagnola_HAP1_v1.0.fa",
         "/Users/sarahcarey/Desktop/hemp_genomes/final_assem/Cannabis_sativa_OttoII_HAP1_v1.0.fa",
         "/Users/sarahcarey/Desktop/hemp_genomes/final_assem/Cannabis_sativa_Futura75_HAP1_v1.0.fa",
         "/Users/sarahcarey/Desktop/hemp_genomes/final_assem/Cannabis_sativa_Uso31_HAP1_v1.0.fa",
         "/Users/sarahcarey/Desktop/hemp_genomes/final_assem/Cannabis_sativa_Carmagnola_HAP2_v1.0.fa",
         "/Users/sarahcarey/Desktop/hemp_genomes/final_assem/Cannabis_sativa_OttoII_HAP2_v1.0.fa",
         "/Users/sarahcarey/Desktop/hemp_genomes/final_assem/Cannabis_sativa_Futura75_HAP2_v1.0.fa",
         "/Users/sarahcarey/Desktop/hemp_genomes/final_assem/Cannabis_sativa_Uso31_HAP2_v1.0.fa")

names(faf) <- c("CsCarH1", "CsOtH1", "CsFuH1", "CsUsH1", 
                "CsCarH2", "CsOtH2", "CsFuH2", "CsUsH2") 

ssl <- sapply(faf, USE.NAMES = T, simplify = F, function(x){
  y <- readDNAStringSet(x)
  return(y)
})

cgtl <- sapply(names(ssl), simplify = F, USE.NAMES = T, function(i){
  print(i)
  find_contigsGapsTelos(
    ssl[[i]],
    teloKmers = c("TTTAGGG"),
    minChrSize = 1e7,
    maxDistBtwTelo = 100,
    minTeloSize = 400,
    minTeloDens = 0.75,
    maxDist2end = 50000)
})
p1 <- plot_contigs(
  cgtList = cgtl,
  palette = viridis::turbo, nColors = 5, nRow = 2)

pdf("Hemp_contigs_telos.pdf", width = 10, height = 8)
p1 <- plot_contigs(
  cgtList = cgtl,
  palette = viridis::turbo, nColors = 5, nRow = 2)
dev.off()


######## Fig. S6, example gene trees #########

#install.packages("ggtree")

library("ggtree")
#v3.14.0

nodeid.tbl_tree <- utils::getFromNamespace("nodeid.tbl_tree", "tidytree")
rootnode.tbl_tree <- utils::getFromNamespace("rootnode.tbl_tree", "tidytree")
offspring.tbl_tree <- utils::getFromNamespace("offspring.tbl_tree", "tidytree")
offspring.tbl_tree_item <- utils::getFromNamespace(".offspring.tbl_tree_item", "tidytree")
child.tbl_tree <- utils::getFromNamespace("child.tbl_tree", "tidytree")
parent.tbl_tree <- utils::getFromNamespace("parent.tbl_tree", "tidytree")

### Shared SDR origin, 4826

tree_file <- read.tree("figure_data/RAxML_bipartitions.OG0004826.tree")

#get node labels
tree_plot <- ggtree(tree_file, color="black", size=1) +
 geom_tiplab(size=2, color="black") +
  geom_nodelab(aes(label=node),size=2, color="red")+
  theme_tree("white")
tree_plot

tree_plot <- ggtree(tree_file, branch.length = "none",
                    size=1, color="gray80") +
  geom_tiplab(size=3.5, color="black") +
  xlim(NA,60) +
  theme_tree("white") +
  geom_nodelab(size=3, color="red")

tree_plot <- collapse(tree_plot, node=91) + 
  geom_point2(aes(subset=(node==91)), shape=23, size=5, fill="black")

tree_plot

ggsave("figure_output/Canna_sharedSDRExample.png", tree_plot, units="in", width=6, height=6, dpi=500,
       device="png")


### Independent captures, 9637

tree_file <- read.tree("figure_data/RAxML_bipartitions.OG0009637.tree")

#get node labels
tree_plot <- ggtree(tree_file, color="black", size=1) +
  geom_tiplab(size=2, color="black") +
  geom_nodelab(aes(label=node),size=2, color="red")+
  theme_tree("white")
tree_plot

tree_plot <- ggtree(tree_file, branch.length = "none",
                    size=1, color="gray80") +
  geom_tiplab(size=3.5, color="black") +
  xlim(NA,60) +
  theme_tree("white") +
  geom_nodelab(size=3, color="red")

tree_plot <- collapse(tree_plot, node=66) + 
  geom_point2(aes(subset=(node==66)), shape=23, size=5, fill="black")

tree_plot <- rotate(tree_plot, node=109)
tree_plot <- rotate(tree_plot, node=120)


tree_plot

ggsave("figure_output/Canna_independentSDRExample.png", tree_plot, units="in", width=6, height=6, dpi=500,
       device="png")



### PAR, 8355


tree_file <- read.tree("figure_data/RAxML_bipartitions.OG0008355.tree")

#get node labels
tree_plot <- ggtree(tree_file, color="black", size=1) +
  geom_tiplab(size=2, color="black") +
  geom_nodelab(aes(label=node),size=2, color="red")+
  theme_tree("white")
tree_plot

tree_plot <- ggtree(tree_file, branch.length = "none",
                    size=1, color="gray80") +
  geom_tiplab(size=3.5, color="black") +
  xlim(NA,60) +
  theme_tree("white") +
  geom_nodelab(size=3, color="red")

tree_plot <- collapse(tree_plot, node=70) + 
  geom_point2(aes(subset=(node==70)), shape=23, size=5, fill="black")

tree_plot

ggsave("figure_output/Canna_PARExample.png", tree_plot, units="in", width=6, height=6, dpi=500,
       device="png")



######## Fig. S8, FT, FD, ACO, and ACS gene trees #########


#get node labels
tree_plot <- ggtree(tree_file, color="black", size=1) +
  geom_tiplab(size=2, color="black") +
  geom_nodelab(aes(label=node),size=2, color="red")+
  theme_tree("white")
tree_plot

### FT, 684

tree_file <- read.tree("figure_data/RAxML_bipartitions.OG0000684.tree")

tree_plot <- ggtree(tree_file, branch.length = "none",
                    size=1, color="gray80") +
  geom_tiplab(size=3.5, color="black") +
  xlim(NA,60) +
  theme_tree("white") +
  geom_nodelab(size=3, color="red")

tree_plot <- collapse(tree_plot, node=459) + 
  geom_point2(aes(subset=(node==459)), shape=23, size=5, fill="black")
tree_plot <- collapse(tree_plot, node=167) + 
  geom_point2(aes(subset=(node==167)), shape=23, size=5, fill="black")

tree_plot

ggsave("figure_output/FT_tree_expanded.png", tree_plot, units="in", width=6, height=10, dpi=500,
       device="png")


### FD, 1983

tree_file <- read.tree("figure_data/RAxML_bipartitions.OG0001913.tree")

#get node labels
tree_plot <- ggtree(tree_file, color="black", size=1) +
  geom_tiplab(size=2, color="black") +
  geom_nodelab(aes(label=node),size=2, color="red")+
  theme_tree("white")
tree_plot


tree_plot <- ggtree(tree_file, branch.length = "none",
                    size=1, color="gray80") +
  geom_tiplab(size=3.5, color="black") +
  xlim(NA,60) +
  theme_tree("white") +
  geom_nodelab(size=3, color="red")

tree_plot <- collapse(tree_plot, node=156) + 
  geom_point2(aes(subset=(node==156)), shape=23, size=5, fill="black")
tree_plot <- collapse(tree_plot, node=168) + 
  geom_point2(aes(subset=(node==168)), shape=23, size=5, fill="black")

tree_plot

ggsave("figure_output/FD_tree.png", tree_plot, units="in", width=6, height=8, dpi=500,
       device="png")


### ACO, 4978

tree_file <- read.tree("figure_data/RAxML_bipartitions.OG0004978.tree")

#get node labels
tree_plot <- ggtree(tree_file, color="black", size=1) +
  geom_tiplab(size=2, color="black") +
  geom_nodelab(aes(label=node),size=2, color="red")+
  theme_tree("white")
tree_plot


tree_plot <- ggtree(tree_file, branch.length = "none",
                    size=1, color="gray80") +
  geom_tiplab(size=3.5, color="black") +
  xlim(NA,60) +
  theme_tree("white") +
  geom_nodelab(size=3, color="red")

tree_plot

ggsave("figure_output/ACO_tree.png", tree_plot, units="in", width=6, height=8, dpi=500,
       device="png")



### ACS, 836

tree_file <- read.tree("figure_data/RAxML_bipartitions.OG0000836.tree")

tree_plot <- ggtree(tree_file, color="black", size=1) +
  geom_tiplab(size=2, color="black") +
  geom_nodelab(aes(label=node),size=2, color="red")+
  theme_tree("white")
tree_plot

tree_plot <- ggtree(tree_file, branch.length = "none",
                    size=1, color="gray80") +
  geom_tiplab(size=3.5, color="black") +
  xlim(NA,60) +
  theme_tree("white") +
  geom_nodelab(size=3, color="red")

tree_plot <- collapse(tree_plot, node=319) + 
  geom_point2(aes(subset=(node==319)), shape=23, size=5, fill="black")

tree_plot

ggsave("figure_output/ACS_tree.png", tree_plot, units="in", width=6, height=8, dpi=500,
       device="png")



################# ANALYSES #################


################# Differential gene expression #################

library("DESeq2")
#v1.46.0

library("ggplot2")
#v3.5.2

countData <- as.matrix(read.csv("figure_data/gene_count_matrix.csv", row.names = "gene_id"))
colData <- read.csv("figure_data/pheno.csv", sep=",", header=TRUE) 


dds <- DESeqDataSetFromMatrix(countData=countData, colData=colData, design=~Tissue_stage_sex)
dds <- DESeq(dds)

# these are the only 3 directly mentioned in the text
res <- results(dds, contrast=c("Tissue_stage_sex","apicalMeristem_L9_female","apicalMeristem_L4_male"))
res2 <- results(dds, contrast=c("Tissue_stage_sex","apicalMeristem_L2_male","apicalMeristem_L4_male"))
res3 <- results(dds, contrast=c("Tissue_stage_sex","leaf_vegetative_female","leaf_flowering_female"))


resOrdered <- res[order(res$pvalue),]
resOrdered2 <- res2[order(res2$pvalue),]
resOrdered3 <- res3[order(res3$pvalue),]


write.csv(as.data.frame(resOrdered), file="results_cannabis_meristem_femaleL9tomaleL4.csv")
write.csv(as.data.frame(resOrdered2), file="results_cannabis_meristem_maleL2tomaleL4.csv")
write.csv(as.data.frame(resOrdered3), file="results_cannabis_leaf_femaleVegetativeToFlowering.csv")


plotMA(res)

vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)
ntd <- normTransform(dds)

plotPCA(vsd, intgroup="Tissue_stage")+ geom_label(aes(label = colData$Sex))




