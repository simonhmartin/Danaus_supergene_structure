
source("~/Research/dev/asynt/asynt.R")

setwd("/home/simon/Research/Danaus_genome/BC_REARRANGEMENTS/alignment_plots")

convert_to_chromosome <- function(alignments, reference_lens, reference_ori, query_lens, query_ori, chrom_name="chrom"){
    a <- alignments
    
    a$Rstart <- get.chrom.pos(scaffold=a$reference, position=a$Rstart,
                              scaf_len=reference_lens, scaf_ori=reference_ori)

    a$Rend <- get.chrom.pos(scaffold=a$reference, position=a$Rend,
                              scaf_len=reference_lens, scaf_ori=reference_ori)
    
    a$Qstart <- get.chrom.pos(scaffold=a$query, position=a$Qstart,
                              scaf_len=query_lens, scaf_ori=query_ori)
    
    a$Qend <- get.chrom.pos(scaffold=a$query, position=a$Qend,
                              scaf_len=query_lens, scaf_ori=query_ori)
    
    a$reference <- chrom_name
    a$query <- chrom_name
    
    a
    }


plex_data <- import.genome(fai_file="../../Dplex_mex/dplex_mex.fa.fai")
plex4_data <- import.genome(fai_file="../../Dplex4/Dapl_Zhan_v3_HiC.chroms.fasta.fai")
MB18102MAT_data <- import.genome(fai_file="../../MB181_trio/MB18102_MAT.fasta.fai")
MB18102PAT_data <- import.genome(fai_file="../../MB181_trio/MB18102_PAT.fasta.fai")
SB211MAT_data <- import.genome(fai_file="../../SB211_trio/SB211_MAT.fasta.fai")
SB211MAT_hifiasm_data <- import.genome(fai_file="../../SB211_trio/SB211_MAT.hifiasm.fasta.fai")
SB211PAT_data <- import.genome(fai_file="../../SB211_trio/SB211_PAT.fasta.fai")
SB211PAT_canu_data <- import.genome(fai_file="../../SB211_trio/SB211_PAT.canu.fasta.fai")
Dchry2.2_data <- import.genome(fai_file="../../Dchry2/Dchry2.2.fa.fai")
Dchry2HAP_data <- import.genome(fai_file="../../Dchry2/Dchry2.haplotigs.fasta.fai")

### regions of interest
### these are regions as defined on the DplexMex scaffold mxdp_6, after rotating it into the correct orientation

regions <- as.data.frame(rbind(c(left=3644267, right=5663000),
                               c(left=5663857, right=7217297),
                               c(left=7217886, right=8499530),
                               c(left=8501152, right=8943270)))



region_cols <- c("#006437", "#a775ee", "#95e79a", "#f7614d")

# region_cols <- c("black","gray30","gray50","gray80")

################################################################################################
###############################  MB18102PAT to plexippus  ######################################
################################################################################################

MB18102PATaln <- import.paf("DplexMex_MB18102PAT.mm2asm20.paf.gz")

reference_contigs <- "mxdp_6"
reference_ori <- c("mxdp_6"="-")

query_contigs <- "tig00000001"
query_ori <- c("tig00000001"="-")

MB18102PATaln <- subset(MB18102PATaln, Rlen>500 & reference == reference_contigs & query %in% query_contigs)

svg("DplexMex_MB18102PAT_chr15_aln.diag.svg", width = 12, height = 12, bg=NA)
par(mfrow = c(1,1), mar = c(2.5,2.5,0,0), xpd=NA)
plot.alignments.diagonal(MB18102PATaln, reference_lens=plex_data$seq_len[reference_contigs], query_lens=MB18102PAT_data$seq_len[query_contigs],
                      reference_ori=reference_ori, query_ori=query_ori, angle_labels=FALSE, labels_cex=1, labels_offset = 3e5,
                      lwd=1, cols=c("#000000","#FF0000"), colour_by ="orientation", xmax=18e6, ymax=18e6, )
dev.off()

svg("DplexMex_MB18102PAT_chr15_aln.svg", width = 12, height = 4, bg=NA)
par(mfrow = c(1,1), mar = c(2,2.5,0,0))
plot.alignments.multi(MB18102PATaln, reference_lens=plex_data$seq_len[reference_contigs], query_lens=MB18102PAT_data$seq_len[query_contigs],
                      reference_ori=reference_ori, query_ori=query_ori, edge_width=0.3,
                      lwd=0.1, cols=c("#000000","#FF0000"), colour_by ="orientation", centre=FALSE, show_labels=TRUE, plot_length=18e6,
                      sigmoid=T, show_outline=FALSE, reference_above=T, gap=0, angle_labels=F, labels_cex=1, labels_offset=0.05)
axis(1, at = seq(0,18e6,1e6), labels=0:18, line=-1)
dev.off()


MB18102PATaln <- convert_to_chromosome(MB18102PATaln,
                                     reference_lens=plex_data$seq_len[reference_contigs], reference_ori=reference_ori,
                                     query_lens=MB18102PAT_data$seq_len[query_contigs], query_ori=query_ori, chrom_name="chr15")

MB18102PATaln <- tidy.alignments(MB18102PATaln)

ref_chrom_len = c("chr15"=sum(plex_data$seq_len[reference_contigs]))
MB18102PAT_chrom_len = c("chr15"=sum(MB18102PAT_data$seq_len[query_contigs]))


#identify synteny blocks
MB18102PATsyn <- get.synteny.blocks.multi(MB18102PATaln, min_subblock_size=200)
MB18102PATsyn <- get.synteny.blocks.multi(MB18102PATsyn, min_subblock_size=2000)
MB18102PATsyn <- get.synteny.blocks.multi(MB18102PATsyn, min_subblock_size=20000)

MB18102PATaln_regions <- list()
MB18102PATsyn_regions <- list()

for (i in 1:nrow(regions)){
    print(i)
    MB18102PATaln_regions[[i]] <- subset(MB18102PATaln, Rstart >= regions[i,"left"] & Rend <= regions[i,"right"])
    
    synblocks <- get.synteny.blocks.multi(MB18102PATaln_regions[[i]], min_subblock_size=200)
    synblocks <- get.synteny.blocks.multi(synblocks, min_subblock_size=2000)
    synblocks <- get.synteny.blocks.multi(synblocks, min_subblock_size=20000)
    
    MB18102PATsyn_regions[[i]] <- synblocks
    
    MB18102PATsyn_regions[[i]]$strand <- ifelse(sign(MB18102PATsyn_regions[[i]]$Rend-MB18102PATsyn_regions[[i]]$Rstart) ==
                                              sign(MB18102PATsyn_regions[[i]]$Qend-MB18102PATsyn_regions[[i]]$Qstart), "+", "-")

    }


# svg("DplexMex_MB18102PAT_chr15_synblocks.wide.svg", width = 9, height = 3, bg=NA)
svg("DplexMex_MB18102PAT_chr15_synblocks.svg", width = 9, height = 4.5, bg=NA)
par(mfrow = c(1,1), mar = c(0,2.5,0,0))
plot.alignments.multi(MB18102PATsyn, reference_lens=ref_chrom_len, query_lens=MB18102PAT_chrom_len,
                      lwd=0.1, cols=c("#000000","#FF0000"), colour_by ="orientation", centre=FALSE, show_labels=FALSE, plot_length=18e6,
                      sigmoid=T, show_outline=FALSE, reference_above=T, gap=0.1)

for (i in 1:nrow(regions)){
    
    draw.horizontal.block.arrow(min(MB18102PATsyn_regions[[i]][,c("Rstart","Rend")]), max(MB18102PATsyn_regions[[i]][,c("Rstart","Rend")]), .05,
                                width=0.05, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
    
    for (j in 1:nrow(MB18102PATsyn_regions[[i]])){
        if (MB18102PATsyn_regions[[i]]$strand[j] == "+"){
            draw.horizontal.block.arrow(min(MB18102PATsyn_regions[[i]][j,c("Qstart","Qend")]), max(MB18102PATsyn_regions[[i]][j,c("Qstart","Qend")]), .95,
                                        width=0.05, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
            }
        else {
            draw.horizontal.block.arrow(max(MB18102PATsyn_regions[[i]][j,c("Qstart","Qend")]), min(MB18102PATsyn_regions[[i]][j,c("Qstart","Qend")]), .95,
                                width=0.05, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
            }
        }
    }

draw.scale.bar(13e6, -0.05, width=1e6, height=0.05, text="1Mb", lwd=1.5)

dev.off()

#output region coordinates
for (i in 1:length(MB18102PATsyn_regions)) MB18102PATsyn_regions[[i]]$Region <- i

write.table(do.call(rbind, MB18102PATsyn_regions), file="MB18102PAT_region_coordinates.tsv",
            sep="\t", quote=F, row.names=F)

################################################################################################
###################################  MB18102MAT to plexippus  ##################################
################################################################################################

MB18102MATaln <- import.paf("DplexMex_MB18102MAT.mm2asm20.paf.gz")

reference_contigs <- "mxdp_6"
reference_ori <- c("mxdp_6"="-")

query_contigs <- c("tig00001608", "tig00001591", "tig00001592")
query_ori <- c("tig00001608"="-", "tig00001591"="+", "tig00001592"="+")

MB18102MATaln <- subset(MB18102MATaln, Rlen>500 & reference == reference_contigs & query %in% query_contigs)

svg("DplexMex_MB18102MAT_chr15_aln.diag.svg", width = 12, height = 12, bg=NA)
par(mfrow = c(1,1), mar = c(2.5,2.5,0,0), xpd=NA)
plot.alignments.diagonal(MB18102MATaln, reference_lens=plex_data$seq_len[reference_contigs], query_lens=MB18102MAT_data$seq_len[query_contigs],
                      reference_ori=reference_ori, query_ori=query_ori, angle_labels=FALSE, labels_cex=1, labels_offset = 3e5,
                      lwd=1, cols=c("#000000","#FF0000"), colour_by ="orientation", xmax=18e6, ymax=18e6, )
dev.off()

svg("DplexMex_MB18102MAT_chr15_aln.svg", width = 12, height = 4, bg=NA)
par(mfrow = c(1,1), mar = c(2,2.5,0,0))
plot.alignments.multi(MB18102MATaln, reference_lens=plex_data$seq_len[reference_contigs], query_lens=MB18102MAT_data$seq_len[query_contigs],
                      reference_ori=reference_ori, query_ori=query_ori, edge_width=0.3,
                      lwd=0.1, cols=c("#000000","#FF0000"), colour_by ="orientation", centre=FALSE, show_labels=TRUE, plot_length=18e6,
                      sigmoid=T, show_outline=FALSE, reference_above=T, gap=0, angle_labels=F, labels_cex=1, labels_offset=0.05)
axis(1, at = seq(0,18e6,1e6), labels=0:18, line=-1)
dev.off()

MB18102MATaln <- convert_to_chromosome(MB18102MATaln,
                                     reference_lens=plex_data$seq_len[reference_contigs], reference_ori=reference_ori,
                                     query_lens=MB18102MAT_data$seq_len[query_contigs], query_ori=query_ori, chrom_name="chr15")

MB18102MATaln <- tidy.alignments(MB18102MATaln)

ref_chrom_len = c("chr15"=sum(plex_data$seq_len[reference_contigs]))
MB18102MAT_chrom_len = c("chr15"=sum(MB18102MAT_data$seq_len[query_contigs]))


#identify synteny blocks
MB18102MATsyn <- get.synteny.blocks.multi(MB18102MATaln, min_subblock_size=200)
MB18102MATsyn <- get.synteny.blocks.multi(MB18102MATsyn, min_subblock_size=2000)
MB18102MATsyn <- get.synteny.blocks.multi(MB18102MATsyn, min_subblock_size=20000)

MB18102MATaln_regions <- list()
MB18102MATsyn_regions <- list()

for (i in 1:nrow(regions)){
    print(i)
    MB18102MATaln_regions[[i]] <- subset(MB18102MATaln, Rstart >= regions[i,"left"] & Rend <= regions[i,"right"])
    
    synblocks <- get.synteny.blocks.multi(MB18102MATaln_regions[[i]], min_subblock_size=200)
    synblocks <- get.synteny.blocks.multi(synblocks, min_subblock_size=2000)
    synblocks <- get.synteny.blocks.multi(synblocks, min_subblock_size=20000)
    
    MB18102MATsyn_regions[[i]] <- synblocks
    
    MB18102MATsyn_regions[[i]]$strand <- ifelse(sign(MB18102MATsyn_regions[[i]]$Rend-MB18102MATsyn_regions[[i]]$Rstart) ==
                                              sign(MB18102MATsyn_regions[[i]]$Qend-MB18102MATsyn_regions[[i]]$Qstart), "+", "-")

    }


# svg("DplexMex_MB18102MAT_chr15_synblocks.wide.svg", width = 9, height = 3, bg=NA)
svg("DplexMex_MB18102MAT_chr15_synblocks.svg", width = 9, height = 4.5, bg=NA)
par(mfrow = c(1,1), mar = c(0,2.5,0,0))
plot.alignments.multi(MB18102MATsyn, reference_lens=ref_chrom_len, query_lens=MB18102MAT_chrom_len,
                      lwd=0.1, cols=c("#000000","#FF0000"), colour_by ="orientation",
                      centre=FALSE, show_labels=FALSE, plot_length=18e6,
                      sigmoid=T, show_outline=FALSE, reference_above=T, gap=0.1)

for (i in 1:nrow(regions)){
    
    draw.horizontal.block.arrow(min(MB18102MATsyn_regions[[i]][,c("Rstart","Rend")]), max(MB18102MATsyn_regions[[i]][,c("Rstart","Rend")]), .05,
                                width=0.05, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
    
    for (j in 1:nrow(MB18102MATsyn_regions[[i]])){
        if (MB18102MATsyn_regions[[i]]$strand[j] == "+"){
            draw.horizontal.block.arrow(min(MB18102MATsyn_regions[[i]][j,c("Qstart","Qend")]), max(MB18102MATsyn_regions[[i]][j,c("Qstart","Qend")]), .95,
                                        width=0.05, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
            }
        else {
            draw.horizontal.block.arrow(max(MB18102MATsyn_regions[[i]][j,c("Qstart","Qend")]), min(MB18102MATsyn_regions[[i]][j,c("Qstart","Qend")]), .95,
                                width=0.05, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
            }
        }
    }


dev.off()

#output region coordinates
for (i in 1:length(MB18102MATsyn_regions)) MB18102MATsyn_regions[[i]]$Region <- i

write.table(do.call(rbind, MB18102MATsyn_regions), file="MB18102MAT_region_coordinates.tsv",
            sep="\t", quote=F, row.names=F)


################################################################################################
######################################  SB211PAT to plexippus  #################################
################################################################################################

SB211PATaln <- import.paf("DplexMex_SB211PAT.mm2asm20.paf.gz")

reference_contigs <- "mxdp_6"
reference_ori <- c("mxdp_6"="-")

query_contigs <- c('h1tg000170l', 'h1tg000359l', 'h1tg000112l', 'h1tg000044l')
query_ori <- c('h1tg000170l'="-", 'h1tg000359l'="+", 'h1tg000112l'="+", 'h1tg000044l'="+")

SB211PATaln <- subset(SB211PATaln, Rlen>500 & reference == reference_contigs & query %in% query_contigs)

svg("DplexMex_SB211PAT_chr15_aln.diag.svg", width = 12, height = 12, bg=NA)
par(mfrow = c(1,1), mar = c(2.5,2.5,0,0), xpd=NA)
plot.alignments.diagonal(SB211PATaln, reference_lens=plex_data$seq_len[reference_contigs], query_lens=SB211PAT_data$seq_len[query_contigs],
                      reference_ori=reference_ori, query_ori=query_ori, angle_labels=TRUE, labels_cex=0.7, labels_offset = 0,
                      lwd=1, cols=c("#000000","#FF0000"), colour_by ="orientation", xmax=18e6, ymax=18e6, )
dev.off()

svg("DplexMex_SB211PAT_chr15_aln.svg", width = 12, height = 4, bg=NA)
par(mfrow = c(1,1), mar = c(2,2.5,0,0))
plot.alignments.multi(SB211PATaln, reference_lens=plex_data$seq_len[reference_contigs], query_lens=SB211PAT_data$seq_len[query_contigs],
                      reference_ori=reference_ori, query_ori=query_ori, edge_width=0.3,
                      lwd=0.1, cols=c("#000000","#FF0000"), colour_by ="orientation", centre=FALSE, show_labels=TRUE, plot_length=18e6,
                      sigmoid=T, show_outline=FALSE, reference_above=T, gap=0, angle_labels=TRUE, labels_cex=0.7, labels_offset=0)
axis(1, at = seq(0,18e6,1e6), labels=0:18, line=-1)
dev.off()

SB211PATaln <- convert_to_chromosome(SB211PATaln,
                                     reference_lens=plex_data$seq_len[reference_contigs], reference_ori=reference_ori,
                                     query_lens=SB211PAT_data$seq_len[query_contigs], query_ori=query_ori, chrom_name="chr15")

ref_chrom_len = c("chr15"=sum(plex_data$seq_len[reference_contigs]))
SB211PAT_chrom_len = c("chr15"=sum(SB211PAT_data$seq_len[query_contigs]))


#identify synteny blocks
SB211PATsyn <- get.synteny.blocks.multi(SB211PATaln, min_subblock_size=200)
SB211PATsyn <- get.synteny.blocks.multi(SB211PATsyn, min_subblock_size=2000)
SB211PATsyn <- get.synteny.blocks.multi(SB211PATsyn, min_subblock_size=20000)

SB211PATaln_regions <- list()
SB211PATsyn_regions <- list()

for (i in 1:nrow(regions)){
    print(i)
    SB211PATaln_regions[[i]] <- subset(SB211PATaln, Rstart >= regions[i,"left"] & Rend <= regions[i,"right"])
    
    synblocks <- get.synteny.blocks.multi(SB211PATaln_regions[[i]], min_subblock_size=200)
    synblocks <- get.synteny.blocks.multi(synblocks, min_subblock_size=2000)
    synblocks <- get.synteny.blocks.multi(synblocks, min_subblock_size=20000)
    
    SB211PATsyn_regions[[i]] <- synblocks
    
    SB211PATsyn_regions[[i]]$strand <- ifelse(sign(SB211PATsyn_regions[[i]]$Rend-SB211PATsyn_regions[[i]]$Rstart) ==
                                              sign(SB211PATsyn_regions[[i]]$Qend-SB211PATsyn_regions[[i]]$Qstart), "+", "-")

    }


# svg("DplexMex_SB211PAT_chr15_synblocks.wide.svg", width = 9, height = 3, bg=NA)
svg("DplexMex_SB211PAT_chr15_synblocks.svg", width = 9, height = 4.5, bg=NA)
par(mfrow = c(1,1), mar = c(0,2.5,0,0))
plot.alignments.multi(SB211PATsyn, reference_lens=ref_chrom_len, query_lens=SB211PAT_chrom_len,
                      lwd=0.1, cols=c("#000000","#FF0000"), colour_by ="orientation", centre=FALSE, show_labels=FALSE, plot_length=18e6,
                      sigmoid=T, show_outline=FALSE, reference_above=T, gap=0.1)

for (i in 1:nrow(regions)){
    
    draw.horizontal.block.arrow(min(SB211PATsyn_regions[[i]][,c("Rstart","Rend")]), max(SB211PATsyn_regions[[i]][,c("Rstart","Rend")]), .05,
                                width=0.05, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
    
    for (j in 1:nrow(SB211PATsyn_regions[[i]])){
        if (SB211PATsyn_regions[[i]]$strand[j] == "+"){
            draw.horizontal.block.arrow(min(SB211PATsyn_regions[[i]][j,c("Qstart","Qend")]), max(SB211PATsyn_regions[[i]][j,c("Qstart","Qend")]), .95,
                                        width=0.05, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
            }
        else {
            draw.horizontal.block.arrow(max(SB211PATsyn_regions[[i]][j,c("Qstart","Qend")]), min(SB211PATsyn_regions[[i]][j,c("Qstart","Qend")]), .95,
                                width=0.05, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
            }
        }
    }


dev.off()

#output region coordinates
for (i in 1:length(SB211PATsyn_regions)) SB211PATsyn_regions[[i]]$Region <- i

write.table(do.call(rbind, SB211PATsyn_regions), file="SB211PAT_region_coordinates.tsv",
            sep="\t", quote=F, row.names=F)



################################################################################################
######################################  SB211MAT to plexippus  #################################
################################################################################################

SB211MATaln <- import.paf("DplexMex_SB211MAT.mm2asm20.paf.gz")

reference_contigs <- "mxdp_6"
reference_ori <- c("mxdp_6"="-")

query_contigs <- c('tig00000001',  'tig00000185')
query_ori <- c('tig00000001'="+",  'tig00000185'="-")

SB211MATaln <- subset(SB211MATaln, Rlen>500 & reference == reference_contigs & query %in% query_contigs)

svg("DplexMex_SB211MAT_chr15_aln.diag.svg", width = 12, height = 12, bg=NA)
par(mfrow = c(1,1), mar = c(2.5,2.5,0,0), xpd=NA)
plot.alignments.diagonal(SB211MATaln, reference_lens=plex_data$seq_len[reference_contigs], query_lens=SB211MAT_data$seq_len[query_contigs],
                      reference_ori=reference_ori, query_ori=query_ori, angle_labels=FALSE, labels_cex=1, labels_offset = 3e5,
                      lwd=1, cols=c("#000000","#FF0000"), colour_by ="orientation", xmax=18e6, ymax=18e6, )
dev.off()

svg("DplexMex_SB211MAT_chr15_aln.svg", width = 12, height = 4, bg=NA)
par(mfrow = c(1,1), mar = c(2,2.5,0,0))
plot.alignments.multi(SB211MATaln, reference_lens=plex_data$seq_len[reference_contigs], query_lens=SB211MAT_data$seq_len[query_contigs],
                      reference_ori=reference_ori, query_ori=query_ori, edge_width=0.3,
                      lwd=0.1, cols=c("#000000","#FF0000"), colour_by ="orientation", centre=FALSE, show_labels=TRUE, plot_length=18e6,
                      sigmoid=T, show_outline=FALSE, reference_above=T, gap=0, angle_labels=F, labels_cex=1, labels_offset=0.05)
axis(1, at = seq(0,18e6,1e6), labels=0:18, line=-1)
dev.off()

SB211MATaln <- convert_to_chromosome(SB211MATaln,
                                     reference_lens=plex_data$seq_len[reference_contigs], reference_ori=reference_ori,
                                     query_lens=SB211MAT_data$seq_len[query_contigs], query_ori=query_ori, chrom_name="chr15")

ref_chrom_len = c("chr15"=sum(plex_data$seq_len[reference_contigs]))
SB211MAT_chrom_len = c("chr15"=sum(SB211MAT_data$seq_len[query_contigs]))


#identify synteny blocks
SB211MATsyn <- get.synteny.blocks.multi(SB211MATaln, min_subblock_size=200)
SB211MATsyn <- get.synteny.blocks.multi(SB211MATsyn, min_subblock_size=2000)
SB211MATsyn <- get.synteny.blocks.multi(SB211MATsyn, min_subblock_size=20000)

SB211MATaln_regions <- list()
SB211MATsyn_regions <- list()

for (i in 1:nrow(regions)){
    print(i)
    SB211MATaln_regions[[i]] <- subset(SB211MATaln, Rstart >= regions[i,"left"] & Rend <= regions[i,"right"])
    
    synblocks <- get.synteny.blocks.multi(SB211MATaln_regions[[i]], min_subblock_size=200)
    synblocks <- get.synteny.blocks.multi(synblocks, min_subblock_size=2000)
    synblocks <- get.synteny.blocks.multi(synblocks, min_subblock_size=20000)
    
    SB211MATsyn_regions[[i]] <- synblocks
    
    SB211MATsyn_regions[[i]]$strand <- ifelse(sign(SB211MATsyn_regions[[i]]$Rend-SB211MATsyn_regions[[i]]$Rstart) ==
                                              sign(SB211MATsyn_regions[[i]]$Qend-SB211MATsyn_regions[[i]]$Qstart), "+", "-")

    }


# svg("DplexMex_SB211MAT_chr15_synblocks.wide.svg", width = 9, height = 3, bg=NA)
svg("DplexMex_SB211MAT_chr15_synblocks.svg", width = 9, height = 4.5, bg=NA)
par(mfrow = c(1,1), mar = c(0,2.5,0,0))
plot.alignments.multi(SB211MATsyn, reference_lens=ref_chrom_len, query_lens=SB211MAT_chrom_len,
                      lwd=0.1, cols=c("#000000","#FF0000"), colour_by ="orientation", centre=FALSE, show_labels=FALSE, plot_length=18e6,
                      sigmoid=T, show_outline=FALSE, reference_above=T, gap=0.1)

for (i in 1:nrow(regions)){
    
    draw.horizontal.block.arrow(min(SB211MATsyn_regions[[i]][,c("Rstart","Rend")]), max(SB211MATsyn_regions[[i]][,c("Rstart","Rend")]), .05,
                                width=0.05, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
    
    for (j in 1:nrow(SB211MATsyn_regions[[i]])){
        if (SB211MATsyn_regions[[i]]$strand[j] == "+"){
            draw.horizontal.block.arrow(min(SB211MATsyn_regions[[i]][j,c("Qstart","Qend")]), max(SB211MATsyn_regions[[i]][j,c("Qstart","Qend")]), .95,
                                        width=0.05, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
            }
        else {
            draw.horizontal.block.arrow(max(SB211MATsyn_regions[[i]][j,c("Qstart","Qend")]), min(SB211MATsyn_regions[[i]][j,c("Qstart","Qend")]), .95,
                                width=0.05, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
            }
        }
    }


dev.off()

#output region coordinates
for (i in 1:length(SB211MATsyn_regions)) SB211MATsyn_regions[[i]]$Region <- i

write.table(do.call(rbind, SB211MATsyn_regions), file="SB211MAT_region_coordinates.tsv",
            sep="\t", quote=F, row.names=F)



################################################################################################
######################################  Dchry2.2 to plexippus  #################################
################################################################################################

Dchry2.2aln <- import.paf("DplexMex_Dchry2.2_minimap2.asm20.paf.gz")

reference_contigs <- "mxdp_6"
reference_ori <- c("mxdp_6"="-")

query_contigs <- 'contig15.1'
query_ori <- c('contig15.1'="+")

Dchry2.2aln <- subset(Dchry2.2aln, Rlen>500 & reference == reference_contigs & query %in% query_contigs)

svg("DplexMex_Dchry2.2_chr15_aln.diag.svg", width = 12, height = 12, bg=NA)
par(mfrow = c(1,1), mar = c(2.5,2.5,0,0), xpd=NA)
plot.alignments.diagonal(Dchry2.2aln, reference_lens=plex_data$seq_len[reference_contigs], query_lens=Dchry2.2_data$seq_len[query_contigs],
                      reference_ori=reference_ori, query_ori=query_ori, angle_labels=FALSE, labels_cex=1, labels_offset = 3e5,
                      lwd=1, cols=c("#000000","#FF0000"), colour_by ="orientation", xmax=18e6, ymax=18e6, )
dev.off()

svg("DplexMex_Dchry2.2_chr15_aln.svg", width = 12, height = 4, bg=NA)
par(mfrow = c(1,1), mar = c(2,2.5,0,0))
plot.alignments.multi(Dchry2.2aln, reference_lens=plex_data$seq_len[reference_contigs], query_lens=Dchry2.2_data$seq_len[query_contigs],
                      reference_ori=reference_ori, query_ori=query_ori, edge_width=0.3,
                      lwd=0.1, cols=c("#000000","#FF0000"), colour_by ="orientation", centre=FALSE, show_labels=TRUE, plot_length=18e6,
                      sigmoid=T, show_outline=FALSE, reference_above=T, gap=0, angle_labels=F, labels_cex=1, labels_offset=0.05)
axis(1, at = seq(0,18e6,1e6), labels=0:18, line=-1)
dev.off()

Dchry2.2aln <- convert_to_chromosome(Dchry2.2aln,
                                     reference_lens=plex_data$seq_len[reference_contigs], reference_ori=reference_ori,
                                     query_lens=Dchry2.2_data$seq_len[query_contigs], query_ori=query_ori, chrom_name="chr15")

ref_chrom_len = c("chr15"=sum(plex_data$seq_len[reference_contigs]))
Dchry2.2_chrom_len = c("chr15"=sum(Dchry2.2_data$seq_len[query_contigs]))


#identify synteny blocks
Dchry2.2syn <- get.synteny.blocks.multi(Dchry2.2aln, min_subblock_size=200)
Dchry2.2syn <- get.synteny.blocks.multi(Dchry2.2syn, min_subblock_size=2000)
Dchry2.2syn <- get.synteny.blocks.multi(Dchry2.2syn, min_subblock_size=20000)

Dchry2.2aln_regions <- list()
Dchry2.2syn_regions <- list()

for (i in 1:nrow(regions)){
    print(i)
    Dchry2.2aln_regions[[i]] <- subset(Dchry2.2aln, Rstart >= regions[i,"left"] & Rend <= regions[i,"right"])
    
    synblocks <- get.synteny.blocks.multi(Dchry2.2aln_regions[[i]], min_subblock_size=200)
    synblocks <- get.synteny.blocks.multi(synblocks, min_subblock_size=2000)
    synblocks <- get.synteny.blocks.multi(synblocks, min_subblock_size=20000)
    
    Dchry2.2syn_regions[[i]] <- synblocks
    
    Dchry2.2syn_regions[[i]]$strand <- ifelse(sign(Dchry2.2syn_regions[[i]]$Rend-Dchry2.2syn_regions[[i]]$Rstart) ==
                                              sign(Dchry2.2syn_regions[[i]]$Qend-Dchry2.2syn_regions[[i]]$Qstart), "+", "-")

    }


# svg("DplexMex_Dchry2.2_chr15_synblocks.wide.svg", width = 9, height = 3, bg=NA)
svg("DplexMex_Dchry2.2_chr15_synblocks.svg", width = 9, height = 4.5, bg=NA)
par(mfrow = c(1,1), mar = c(0,2.5,0,0))
plot.alignments.multi(Dchry2.2syn, reference_lens=ref_chrom_len, query_lens=Dchry2.2_chrom_len,
                      lwd=0.1, cols=c("#000000","#FF0000"), colour_by ="orientation", centre=FALSE, show_labels=FALSE, plot_length=18e6,
                      sigmoid=T, show_outline=FALSE, reference_above=T, gap=0.1)

for (i in 1:nrow(regions)){
    
    draw.horizontal.block.arrow(min(Dchry2.2syn_regions[[i]][,c("Rstart","Rend")]), max(Dchry2.2syn_regions[[i]][,c("Rstart","Rend")]), .05,
                                width=0.05, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
    
    for (j in 1:nrow(Dchry2.2syn_regions[[i]])){
        if (Dchry2.2syn_regions[[i]]$strand[j] == "+"){
            draw.horizontal.block.arrow(min(Dchry2.2syn_regions[[i]][j,c("Qstart","Qend")]), max(Dchry2.2syn_regions[[i]][j,c("Qstart","Qend")]), .95,
                                        width=0.05, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
            }
        else {
            draw.horizontal.block.arrow(max(Dchry2.2syn_regions[[i]][j,c("Qstart","Qend")]), min(Dchry2.2syn_regions[[i]][j,c("Qstart","Qend")]), .95,
                                width=0.05, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
            }
        }
    }

draw.scale.bar(13e6, -0.05, width=1e6, height=0.05, text="1Mb", lwd=1.5)

dev.off()

#output region coordinates
for (i in 1:length(Dchry2.2syn_regions)) Dchry2.2syn_regions[[i]]$Region <- i

write.table(do.call(rbind, Dchry2.2syn_regions), file="Dchry2.2_region_coordinates.tsv",
            sep="\t", quote=F, row.names=F)



################################################################################################
###################################  Dchry2HAP to plexippus  ################################
################################################################################################

Dchry2HAPaln <- import.paf("DplexMex_Dchry2hap.mm2asm20.paf.gz")

reference_contigs <- "mxdp_6"
reference_ori <- c("mxdp_6"="-")

query_contigs <- c("000397F")
query_ori <- c("000397F"="-")

Dchry2HAPaln <- subset(Dchry2HAPaln, Rlen>500 & reference == reference_contigs & query %in% query_contigs)

svg("DplexMex_Dchry2HAP_chr15_aln.diag.svg", width = 12, height = 12, bg=NA)
par(mfrow = c(1,1), mar = c(2.5,2.5,0,0), xpd=NA)
plot.alignments.diagonal(Dchry2HAPaln, reference_lens=plex_data$seq_len[reference_contigs], query_lens=Dchry2HAP_data$seq_len[query_contigs],
                      reference_ori=reference_ori, query_ori=query_ori, angle_labels=FALSE, labels_cex=1, labels_offset = 3e5,
                      lwd=1, cols=c("#000000","#FF0000"), colour_by ="orientation", xmax=18e6, ymax=18e6, )
dev.off()

svg("DplexMex_Dchry2HAP_chr15_aln.svg", width = 12, height = 4, bg=NA)
par(mfrow = c(1,1), mar = c(2,2.5,0,0))
plot.alignments.multi(Dchry2HAPaln, reference_lens=plex_data$seq_len[reference_contigs], query_lens=Dchry2HAP_data$seq_len[query_contigs],
                      reference_ori=reference_ori, query_ori=query_ori, edge_width=0.3,
                      lwd=0.1, cols=c("#000000","#FF0000"), colour_by ="orientation", centre=FALSE, show_labels=TRUE, plot_length=18e6,
                      sigmoid=T, show_outline=FALSE, reference_above=T, gap=0, angle_labels=F, labels_cex=1, labels_offset=0.05)
axis(1, at = seq(0,18e6,1e6), labels=0:18, line=-1)
dev.off()


################################################################################################
#################################  MB18102MAT to MB18102PAT  ###################################
################################################################################################

MB18102MATtoPATaln <- import.paf("MB18102PAT_MB18102MAT_minimap2.asm10.paf.gz", additional_fields="dv")

reference_contigs <- "tig00000001"
reference_ori <- c("tig00000001"="-")

query_contigs <- c("tig00001608", "tig00001591", "tig00001592")
query_ori <- c("tig00001608"="-", "tig00001591"="+", "tig00001592"="+")


MB18102MATtoPATaln <- subset(MB18102MATtoPATaln, reference %in% reference_contigs & query %in% query_contigs)

MB18102MATtoPATaln <- subset(MB18102MATtoPATaln, dv <= 0.2 & Rlen > 2000)

svg("MB18102PAT_MB18102MAT_chr15.svg", width = 15, height = 5)
par(mfrow = c(1,1), mar = c(0,0,0,0))
plot.alignments.multi(MB18102MATtoPATaln, reference_lens=MB18102PAT_data$seq_len, query_lens=MB18102MAT_data$seq_len, reference_ori=reference_ori, query_ori=query_ori,
                      lwd=0.1, cols=c("#000000","#ff0000"), colour_by ="orientation", centre=FALSE, show_labels=TRUE, plot_length=18e6,
                      sigmoid=T, show_outline=FALSE, reference_above=T)
dev.off()

#as chromosomes
MB18102MATtoPATaln <- convert_to_chromosome(MB18102MATtoPATaln,
                                     reference_lens=MB18102PAT_data$seq_len[reference_contigs], reference_ori=reference_ori,
                                     query_lens=MB18102MAT_data$seq_len[query_contigs], query_ori=query_ori, chrom_name="chr15")


#identify synteny blocks
MB18102MATtoPATsyn <- get.synteny.blocks.multi(MB18102MATtoPATaln, min_subblock_size=200, min_block_size=50000)
MB18102MATtoPATsyn <- get.synteny.blocks.multi(MB18102MATtoPATsyn, min_subblock_size=2000, min_block_size=50000)
MB18102MATtoPATsyn <- get.synteny.blocks.multi(MB18102MATtoPATsyn, min_subblock_size=20000, min_block_size=50000)


svg("MB18102PAT_MB18102MAT_chr15_synblocks.svg", width = 7, height = 3.5, bg=NA)
par(mfrow = c(1,1), mar = c(0,2.5,0,0))
plot.alignments.multi(MB18102MATtoPATsyn, reference_lens=MB18102PAT_chrom_len, query_lens=MB18102MAT_chrom_len,
                      lwd=0.1, cols=c("#000000","#ff0000"), colour_by ="orientation", centre=FALSE, show_labels=FALSE, plot_length=18e6,
                      sigmoid=T, show_outline=FALSE, reference_above=T, gap=0.1)



for (i in 1:nrow(regions)){
    
    for (j in 1:nrow(MB18102PATsyn_regions[[i]])){
        if (MB18102PATsyn_regions[[i]]$strand[j] == "+"){
            draw.horizontal.block.arrow(min(MB18102PATsyn_regions[[i]][j,c("Qstart","Qend")]), max(MB18102PATsyn_regions[[i]][j,c("Qstart","Qend")]), .05,
                                        width=0.05, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
            }
        else {
            draw.horizontal.block.arrow(max(MB18102PATsyn_regions[[i]][j,c("Qstart","Qend")]), min(MB18102PATsyn_regions[[i]][j,c("Qstart","Qend")]), .05,
                                width=0.05, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
            }
        }
    
    for (j in 1:nrow(MB18102MATsyn_regions[[i]])){
        if (MB18102MATsyn_regions[[i]]$strand[j] == "+"){
            draw.horizontal.block.arrow(min(MB18102MATsyn_regions[[i]][j,c("Qstart","Qend")]), max(MB18102MATsyn_regions[[i]][j,c("Qstart","Qend")]), .95,
                                        width=0.05, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
            }
        else {
            draw.horizontal.block.arrow(max(MB18102MATsyn_regions[[i]][j,c("Qstart","Qend")]), min(MB18102MATsyn_regions[[i]][j,c("Qstart","Qend")]), .95,
                                width=0.05, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
            }
        }
    }



dev.off()



################################################################################################
###################### SB211PAT (orientis) to MB18102MAT (chrysippus)  #########################
################################################################################################


SB211PATtoMB18102MATaln <- import.paf("MB18102MAT_SB211PAT.mm2asm10.paf.gz", additional_fields="dv")

reference_contigs <- c("tig00001608", "tig00001591", "tig00001592")
reference_ori <- c("tig00001608"="-", "tig00001591"="+", "tig00001592"="+")

query_contigs <- c('h1tg000170l', 'h1tg000359l', 'h1tg000112l', 'h1tg000044l')
query_ori <- c('h1tg000170l'="-", 'h1tg000359l'="+", 'h1tg000112l'="+", 'h1tg000044l'="+")

SB211PATtoMB18102MATaln <- subset(SB211PATtoMB18102MATaln, reference %in% reference_contigs & query %in% query_contigs)

SB211PATtoMB18102MATaln <- subset(SB211PATtoMB18102MATaln, dv <= 0.2 & Rlen > 2000)


svg("MB18102MAT_SB211PAT_chr15.svg", width = 15, height = 5)
par(mfrow = c(1,1), mar = c(0,0,0,0))
plot.alignments.multi(SB211PATtoMB18102MATaln, reference_lens=MB18102MAT_data$seq_len[reference_contigs], query_lens=SB211PAT_data$seq_len[query_contigs],
                      reference_ori=reference_ori, query_ori=query_ori, no_reorder=TRUE,
                      lwd=0.1, cols=c("#000000","#ff0000"), colour_by ="orientation", centre=FALSE, show_labels=TRUE, plot_length=18e6,
                      sigmoid=T, show_outline=FALSE, reference_above=T)
dev.off()

#as chromosomes
SB211PATtoMB18102MATaln <- convert_to_chromosome(SB211PATtoMB18102MATaln,
                                     reference_lens=MB18102MAT_data$seq_len[reference_contigs], reference_ori=reference_ori,
                                     query_lens=SB211PAT_data$seq_len[query_contigs], query_ori=query_ori, chrom_name="chr15")


#identify synteny blocks
SB211PATtoMB18102MATsyn <- get.synteny.blocks.multi(SB211PATtoMB18102MATaln, min_subblock_size=200, min_block_size=50000)
SB211PATtoMB18102MATsyn <- get.synteny.blocks.multi(SB211PATtoMB18102MATsyn, min_subblock_size=2000, min_block_size=50000)
SB211PATtoMB18102MATsyn <- get.synteny.blocks.multi(SB211PATtoMB18102MATsyn, min_subblock_size=20000, min_block_size=50000)


svg("MB18102MAT_SB211PAT_chr15_synblocks.svg", width = 7, height = 3.5, bg=NA)
par(mfrow = c(1,1), mar = c(0,2.5,0,0))
plot.alignments.multi(SB211PATtoMB18102MATsyn, reference_lens=MB18102MAT_chrom_len, query_lens=SB211PAT_chrom_len,
                      lwd=0.1, cols=c("#000000","#ff0000"), colour_by ="orientation", centre=FALSE, show_labels=FALSE, plot_length=18e6,
                      sigmoid=T, show_outline=FALSE, reference_above=T, gap=0.1)


for (i in 1:nrow(regions)){
    
    for (j in 1:nrow(MB18102MATsyn_regions[[i]])){
        if (MB18102MATsyn_regions[[i]]$strand[j] == "+"){
            draw.horizontal.block.arrow(min(MB18102MATsyn_regions[[i]][j,c("Qstart","Qend")]), max(MB18102MATsyn_regions[[i]][j,c("Qstart","Qend")]), .05,
                                        width=0.05, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
            }
        else {
            draw.horizontal.block.arrow(max(MB18102MATsyn_regions[[i]][j,c("Qstart","Qend")]), min(MB18102MATsyn_regions[[i]][j,c("Qstart","Qend")]), .05,
                                width=0.05, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
            }
        }
    
    for (j in 1:nrow(SB211PATsyn_regions[[i]])){
        if (SB211PATsyn_regions[[i]]$strand[j] == "+"){
            draw.horizontal.block.arrow(min(SB211PATsyn_regions[[i]][j,c("Qstart","Qend")]), max(SB211PATsyn_regions[[i]][j,c("Qstart","Qend")]), .95,
                                        width=0.05, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
            }
        else {
            draw.horizontal.block.arrow(max(SB211PATsyn_regions[[i]][j,c("Qstart","Qend")]), min(SB211PATsyn_regions[[i]][j,c("Qstart","Qend")]), .95,
                                width=0.05, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
            }
        }
    }



dev.off()


################################################################################################
#################################  MB18102PAT to SB211PAT  ###################################
################################################################################################

MB18102PATtoSB211PATaln <- import.paf("SB211PAT_MB18102PAT_minimap2.asm10.paf.gz", additional_fields="dv")

reference_contigs <- c('h1tg000170l', 'h1tg000359l', 'h1tg000112l', 'h1tg000044l')
reference_ori <- c('h1tg000170l'="-", 'h1tg000359l'="+", 'h1tg000112l'="+", 'h1tg000044l'="+")

query_contigs <- "tig00000001"
query_ori <- c("tig00000001"="-")

MB18102PATtoSB211PATaln <- subset(MB18102PATtoSB211PATaln, reference %in% reference_contigs & query %in% query_contigs)

MB18102PATtoSB211PATaln <- subset(MB18102PATtoSB211PATaln, dv <= 0.2 & Rlen > 2000)

svg("SB211PAT_MB18102PAT_chr15.svg", width = 15, height = 5)
par(mfrow = c(1,1), mar = c(0,0,0,0))
plot.alignments.multi(MB18102PATtoSB211PATaln, reference_lens=SB211PAT_data$seq_len[reference_contigs], query_lens=MB18102PAT_data$seq_len,
                      reference_ori=reference_ori, query_ori=query_ori,
                      lwd=0.1, cols=c("#000000","#ff0000"), colour_by ="orientation", centre=FALSE, show_labels=TRUE, plot_length=18e6,
                      sigmoid=T, show_outline=FALSE, reference_above=T)
dev.off()


#as chromosomes
MB18102PATtoSB211PATaln <- convert_to_chromosome(MB18102PATtoSB211PATaln,
                                     reference_lens=SB211PAT_data$seq_len[reference_contigs], reference_ori=reference_ori,
                                     query_lens=MB18102PAT_data$seq_len[query_contigs], query_ori=query_ori, chrom_name="chr15")


#identify synteny blocks
MB18102PATtoSB211PATsyn <- get.synteny.blocks.multi(MB18102PATtoSB211PATaln, min_subblock_size=200, min_block_size=50000)
MB18102PATtoSB211PATsyn <- get.synteny.blocks.multi(MB18102PATtoSB211PATsyn, min_subblock_size=2000, min_block_size=50000)
MB18102PATtoSB211PATsyn <- get.synteny.blocks.multi(MB18102PATtoSB211PATsyn, min_subblock_size=20000, min_block_size=50000)


svg("SB211PAT_MB18102PAT_chr15_synblocks.svg", width = 7, height = 3.5, bg=NA)
par(mfrow = c(1,1), mar = c(0,2.5,0,0))
plot.alignments.multi(MB18102PATtoSB211PATsyn, reference_lens=SB211PAT_chrom_len, query_lens=MB18102PAT_chrom_len,
                      lwd=0.1, cols=c("#000000","#ff0000"), colour_by ="orientation", centre=FALSE, show_labels=FALSE, plot_length=18e6,
                      sigmoid=T, show_outline=FALSE, reference_above=T, gap=0.1)



for (i in 1:nrow(regions)){
    
    for (j in 1:nrow(MB18102PATsyn_regions[[i]])){
        if (MB18102PATsyn_regions[[i]]$strand[j] == "+"){
            draw.horizontal.block.arrow(min(MB18102PATsyn_regions[[i]][j,c("Qstart","Qend")]), max(MB18102PATsyn_regions[[i]][j,c("Qstart","Qend")]), .95,
                                        width=0.05, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
            }
        else {
            draw.horizontal.block.arrow(max(MB18102PATsyn_regions[[i]][j,c("Qstart","Qend")]), min(MB18102PATsyn_regions[[i]][j,c("Qstart","Qend")]), .95,
                                width=0.05, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
            }
        }
    
    for (j in 1:nrow(SB211PATsyn_regions[[i]])){
        if (SB211PATsyn_regions[[i]]$strand[j] == "+"){
            draw.horizontal.block.arrow(min(SB211PATsyn_regions[[i]][j,c("Qstart","Qend")]), max(SB211PATsyn_regions[[i]][j,c("Qstart","Qend")]), .05,
                                        width=0.05, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
            }
        else {
            draw.horizontal.block.arrow(max(SB211PATsyn_regions[[i]][j,c("Qstart","Qend")]), min(SB211PATsyn_regions[[i]][j,c("Qstart","Qend")]), .05,
                                width=0.05, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
            }
        }
    }


dev.off()



################################################################################################
###################################  MB18102MAT to Dchry2.2  ###################################
################################################################################################

MB18102MATtoDchry2.2aln <- import.paf("Dchry2.2_MB18102MAT_mm2asm10.paf.gz", additional_fields="dv")

reference_contigs <- "contig15.1"
reference_ori <- c("contig15.1"="+")

query_contigs <- c("tig00001608", "tig00001591", "tig00001592")
query_ori <- c("tig00001608"="-", "tig00001591"="+", "tig00001592"="+")


MB18102MATtoDchry2.2aln <- subset(MB18102MATtoDchry2.2aln, reference %in% reference_contigs & query %in% query_contigs)

MB18102MATtoDchry2.2aln <- subset(MB18102MATtoDchry2.2aln, dv <= 0.2 & Rlen > 2000)

svg("Dchry2.2_MB18102MAT_chr15.svg", width = 15, height = 5)
par(mfrow = c(1,1), mar = c(0,0,0,0))
plot.alignments.multi(MB18102MATtoDchry2.2aln, reference_lens=Dchry2.2_data$seq_len, query_lens=MB18102MAT_data$seq_len, reference_ori=reference_ori, query_ori=query_ori,
                      lwd=0.1, cols=c("#000000","#ff0000"), colour_by ="orientation", centre=FALSE, show_labels=TRUE, plot_length=18e6,
                      sigmoid=T, show_outline=FALSE, reference_above=T)
dev.off()

#as chromosomes
MB18102MATtoDchry2.2aln <- convert_to_chromosome(MB18102MATtoDchry2.2aln,
                                     reference_lens=Dchry2.2_data$seq_len[reference_contigs], reference_ori=reference_ori,
                                     query_lens=MB18102MAT_data$seq_len[query_contigs], query_ori=query_ori, chrom_name="chr15")


#identify synteny blocks
MB18102MATtoDchry2.2syn <- get.synteny.blocks.multi(MB18102MATtoDchry2.2aln, min_subblock_size=200, min_block_size=50000)
MB18102MATtoDchry2.2syn <- get.synteny.blocks.multi(MB18102MATtoDchry2.2syn, min_subblock_size=2000, min_block_size=50000)
MB18102MATtoDchry2.2syn <- get.synteny.blocks.multi(MB18102MATtoDchry2.2syn, min_subblock_size=20000, min_block_size=50000)


svg("Dchry2.2_MB18102MAT_chr15_synblocks.svg", width = 7, height = 3.5, bg=NA)
par(mfrow = c(1,1), mar = c(0,2.5,0,0))
plot.alignments.multi(MB18102MATtoDchry2.2syn, reference_lens=Dchry2.2_chrom_len, query_lens=MB18102MAT_chrom_len,
                      lwd=0.1, cols=c("#000000","#ff0000"), colour_by ="orientation", centre=FALSE, show_labels=FALSE, plot_length=18e6,
                      sigmoid=T, show_outline=FALSE, reference_above=T, gap=0.1)



for (i in 1:nrow(regions)){
    
    for (j in 1:nrow(Dchry2.2syn_regions[[i]])){
        if (Dchry2.2syn_regions[[i]]$strand[j] == "+"){
            draw.horizontal.block.arrow(min(Dchry2.2syn_regions[[i]][j,c("Qstart","Qend")]), max(Dchry2.2syn_regions[[i]][j,c("Qstart","Qend")]), .05,
                                        width=0.05, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
            }
        else {
            draw.horizontal.block.arrow(max(Dchry2.2syn_regions[[i]][j,c("Qstart","Qend")]), min(Dchry2.2syn_regions[[i]][j,c("Qstart","Qend")]), .05,
                                width=0.05, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
            }
        }
    
    for (j in 1:nrow(MB18102MATsyn_regions[[i]])){
        if (MB18102MATsyn_regions[[i]]$strand[j] == "+"){
            draw.horizontal.block.arrow(min(MB18102MATsyn_regions[[i]][j,c("Qstart","Qend")]), max(MB18102MATsyn_regions[[i]][j,c("Qstart","Qend")]), .95,
                                        width=0.05, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
            }
        else {
            draw.horizontal.block.arrow(max(MB18102MATsyn_regions[[i]][j,c("Qstart","Qend")]), min(MB18102MATsyn_regions[[i]][j,c("Qstart","Qend")]), .95,
                                width=0.05, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
            }
        }
    }



dev.off()


################################################################################################
#################################  Dchry2.2 to SB211PAT  ###################################
################################################################################################

Dchry2.2toSB211PATaln <- import.paf("SB211PAT_Dchry2.2_mm2asm10.paf.gz", additional_fields="dv")

reference_contigs <- c('h1tg000170l', 'h1tg000359l', 'h1tg000112l', 'h1tg000044l')
reference_ori <- c('h1tg000170l'="-", 'h1tg000359l'="+", 'h1tg000112l'="+", 'h1tg000044l'="+")

query_contigs <- "contig15.1"
query_ori <- c("contig15.1"="+")

Dchry2.2toSB211PATaln <- subset(Dchry2.2toSB211PATaln, reference %in% reference_contigs & query %in% query_contigs)

Dchry2.2toSB211PATaln <- subset(Dchry2.2toSB211PATaln, dv <= 0.2 & Rlen > 2000)

svg("SB211PAT_Dchry2.2_chr15.svg", width = 15, height = 5)
par(mfrow = c(1,1), mar = c(0,0,0,0))
plot.alignments.multi(Dchry2.2toSB211PATaln, reference_lens=SB211PAT_data$seq_len[reference_contigs], query_lens=Dchry2.2_data$seq_len,
                      reference_ori=reference_ori, query_ori=query_ori,
                      lwd=0.1, cols=c("#000000","#ff0000"), colour_by ="orientation", centre=FALSE, show_labels=TRUE, plot_length=18e6,
                      sigmoid=T, show_outline=FALSE, reference_above=T)
dev.off()


#as chromosomes
Dchry2.2toSB211PATaln <- convert_to_chromosome(Dchry2.2toSB211PATaln,
                                     reference_lens=SB211PAT_data$seq_len[reference_contigs], reference_ori=reference_ori,
                                     query_lens=Dchry2.2_data$seq_len[query_contigs], query_ori=query_ori, chrom_name="chr15")


#identify synteny blocks
Dchry2.2toSB211PATsyn <- get.synteny.blocks.multi(Dchry2.2toSB211PATaln, min_subblock_size=200, min_block_size=50000)
Dchry2.2toSB211PATsyn <- get.synteny.blocks.multi(Dchry2.2toSB211PATsyn, min_subblock_size=2000, min_block_size=50000)
Dchry2.2toSB211PATsyn <- get.synteny.blocks.multi(Dchry2.2toSB211PATsyn, min_subblock_size=20000, min_block_size=50000)


svg("SB211PAT_Dchry2.2_chr15_synblocks.svg", width = 7, height = 3.5, bg=NA)
par(mfrow = c(1,1), mar = c(0,2.5,0,0))
plot.alignments.multi(Dchry2.2toSB211PATsyn, reference_lens=SB211PAT_chrom_len, query_lens=Dchry2.2_chrom_len,
                      lwd=0.1, cols=c("#000000","#ff0000"), colour_by ="orientation", centre=FALSE, show_labels=FALSE, plot_length=18e6,
                      sigmoid=T, show_outline=FALSE, reference_above=T, gap=0.1)



for (i in 1:nrow(regions)){
    
    for (j in 1:nrow(Dchry2.2syn_regions[[i]])){
        if (Dchry2.2syn_regions[[i]]$strand[j] == "+"){
            draw.horizontal.block.arrow(min(Dchry2.2syn_regions[[i]][j,c("Qstart","Qend")]), max(Dchry2.2syn_regions[[i]][j,c("Qstart","Qend")]), .95,
                                        width=0.05, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
            }
        else {
            draw.horizontal.block.arrow(max(Dchry2.2syn_regions[[i]][j,c("Qstart","Qend")]), min(Dchry2.2syn_regions[[i]][j,c("Qstart","Qend")]), .95,
                                width=0.05, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
            }
        }
    
    for (j in 1:nrow(SB211PATsyn_regions[[i]])){
        if (SB211PATsyn_regions[[i]]$strand[j] == "+"){
            draw.horizontal.block.arrow(min(SB211PATsyn_regions[[i]][j,c("Qstart","Qend")]), max(SB211PATsyn_regions[[i]][j,c("Qstart","Qend")]), .05,
                                        width=0.05, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
            }
        else {
            draw.horizontal.block.arrow(max(SB211PATsyn_regions[[i]][j,c("Qstart","Qend")]), min(SB211PATsyn_regions[[i]][j,c("Qstart","Qend")]), .05,
                                width=0.05, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
            }
        }
    }


dev.off()

################################################################################################
##############################  MB18102PAT to Dchry2.2  ###################################
################################################################################################

MB18102PATtoDchry2.2aln <- import.paf("Dchry2.2_MB18102PAT_mm2asm10.paf.gz", additional_fields="dv")

reference_contigs <- "contig15.1"
reference_ori <- c("contig15.1"="+")

query_contigs <- "tig00000001"
query_ori <- c("tig00000001"="-")


MB18102PATtoDchry2.2aln <- subset(MB18102PATtoDchry2.2aln, reference %in% reference_contigs & query %in% query_contigs)

MB18102PATtoDchry2.2aln <- subset(MB18102PATtoDchry2.2aln, dv <= 0.2 & Rlen > 2000)

svg("Dchry2.2_MB18102PAT_chr15_aln.diag.svg", width = 6, height = 6, bg=NA)
par(mfrow = c(1,1), mar = c(2.5,2.5,0,0), xpd=NA)
plot.alignments.diagonal(MB18102PATtoDchry2.2aln, reference_lens=Dchry2.2_data$seq_len[reference_contigs], query_lens=MB18102PAT_data$seq_len[query_contigs],
                      reference_ori=reference_ori, query_ori=query_ori, labels_cex=1, angle_labels=F, labels_offset=3e5,
                      lwd=1, cols=c("#000000","#FF0000"), colour_by ="orientation", xmax=18e6, ymax=18e6)
dev.off()

svg("Dchry2.2_MB18102PAT_chr15_aln.svg", width = 12, height = 4, bg=NA)
par(mfrow = c(1,1), mar = c(2,2.5,0,0))
plot.alignments.multi(MB18102PATtoDchry2.2aln, reference_lens=Dchry2.2_data$seq_len[reference_contigs], query_lens=MB18102PAT_data$seq_len[query_contigs],
                      reference_ori=reference_ori, query_ori=query_ori, edge_width=0.5,
                      lwd=0.1, cols=c("#000000","#FF0000"), colour_by ="orientation", centre=FALSE, show_labels=TRUE, plot_length=18e6,
                      sigmoid=T, show_outline=FALSE, reference_above=T, gap=0)
axis(1, at = seq(0,18e6,1e6), labels=0:18, line=-1)
dev.off()


################################################################################################
##############################  MB18102MAT to Dchry2HAP  ###################################
################################################################################################

MB18102MATtoDchry2HAPaln <- import.paf("Dchry2HAP_MB18102MAT_mm2asm10.paf.gz", additional_fields="dv")

reference_contigs <- "000397F"
reference_ori <- c("000397F"="-")

query_contigs <- c("tig00001608", "tig00001591", "tig00001592")
query_ori <- c("tig00001608"="-", "tig00001591"="+", "tig00001592"="+")


MB18102MATtoDchry2HAPaln <- subset(MB18102MATtoDchry2HAPaln, reference %in% reference_contigs & query %in% query_contigs)

MB18102MATtoDchry2HAPaln <- subset(MB18102MATtoDchry2HAPaln, dv <= 0.2 & Rlen > 2000)

svg("Dchry2HAP_MB18102MAT_chr15_aln.diag.svg", width = 6, height = 6, bg=NA)
par(mfrow = c(1,1), mar = c(2.5,2.5,0,0), xpd=NA)
plot.alignments.diagonal(MB18102MATtoDchry2HAPaln, reference_lens=Dchry2HAP_data$seq_len[reference_contigs], query_lens=MB18102MAT_data$seq_len[query_contigs],
                      reference_ori=reference_ori, query_ori=query_ori, labels_cex=1, angle_labels=F, labels_offset=3e5,
                      lwd=1, cols=c("#000000","#FF0000"), colour_by ="orientation", xmax=18e6, ymax=18e6)
dev.off()

svg("Dchry2HAP_MB18102MAT_chr15_aln.svg", width = 12, height = 4, bg=NA)
par(mfrow = c(1,1), mar = c(2,2.5,0,0))
plot.alignments.multi(MB18102MATtoDchry2HAPaln, reference_lens=Dchry2HAP_data$seq_len[reference_contigs], query_lens=MB18102MAT_data$seq_len[query_contigs],
                      reference_ori=reference_ori, query_ori=query_ori, edge_width=0.5,
                      lwd=0.1, cols=c("#000000","#FF0000"), colour_by ="orientation", centre=FALSE, show_labels=TRUE, plot_length=18e6,
                      sigmoid=T, show_outline=FALSE, reference_above=T, gap=0)
axis(1, at = seq(0,18e6,1e6), labels=0:18, line=-1)
dev.off()

################################################################################################
##############################  SB211MAT hifiasm vs canu  ###################################
################################################################################################

SB211MAT_hifiasm_canu_aln <- import.paf("SB211MAT_hifiasm_vs_canu.mm2asm10.paf.gz", additional_fields="dv")

reference_contigs <- c('h1tg000114l',"h1tg000429l",    "h1tg000018l",    "h1tg000280l",    "h1tg000474l",    "h1tg000041l",    "h1tg000123l",    "h1tg000587l")
reference_ori <- c('h1tg000114l'="+","h1tg000429l"="-","h1tg000018l"="+","h1tg000280l"="+","h1tg000474l"="-","h1tg000041l"="+","h1tg000123l"="-","h1tg000587l"="-")

query_contigs <- c('tig00000001',  'tig00000185')
query_ori <- c('tig00000001'="+",  'tig00000185'="+")

SB211MAT_hifiasm_canu_aln <- subset(SB211MAT_hifiasm_canu_aln, reference %in% reference_contigs & query %in% query_contigs)

SB211MAT_hifiasm_canu_aln <- subset(SB211MAT_hifiasm_canu_aln, dv <= 0.2 & Rlen > 2000)

svg("SB211MAT_hifiasm_canu_aln.diag.svg", width = 6, height = 6, bg=NA)
par(mfrow = c(1,1), mar = c(2.5,2.5,0,0), xpd=NA)
plot.alignments.diagonal(SB211MAT_hifiasm_canu_aln, reference_lens=SB211MAT_hifiasm_data$seq_len[reference_contigs], query_lens=SB211MAT_data$seq_len[query_contigs],
                      reference_ori=reference_ori, query_ori=query_ori, labels_cex=0.8, angle_labels=T, labels_offset=3e5,
                      lwd=1, cols=c("#000000","#FF0000"), colour_by ="orientation", xmax=18e6, ymax=18e6)
dev.off()


################################################################################################
##############################  SB211MAT hifiasm vs canu  ###################################
################################################################################################

SB211PAT_hifiasm_canu_aln <- import.paf("SB211PAT_hifiasm_vs_canu.mm2asm10.paf.gz", additional_fields="dv")

reference_contigs <- c('h1tg000170l', 'h1tg000359l', 'h1tg000112l', 'h1tg000044l')
reference_ori <- c('h1tg000170l'="-", 'h1tg000359l'="+", 'h1tg000112l'="+", 'h1tg000044l'="+")

query_contigs <- c("tig00004103", "tig00000150", "tig00000176", "tig00004106", "tig00000336",
                   "tig00000540", "tig00000617", "tig00000632", "tig00000658", "tig00001154",
                   "tig00001001", "tig00001268", "tig00001468", "tig00002919")

query_ori <- c(tig00004103="+", tig00000150="-", tig00000176="-", tig00004106="+",
                tig00000336="+", tig00000540="-", tig00000617="-", tig00000632="-",
                tig00000658="-", tig00001154="+", tig00001001="-", tig00001268="+",
                tig00001468="+", tig00002919="-")

SB211PAT_hifiasm_canu_aln <- subset(SB211PAT_hifiasm_canu_aln, reference %in% reference_contigs & query %in% query_contigs)

SB211PAT_hifiasm_canu_aln <- subset(SB211PAT_hifiasm_canu_aln, dv <= 0.2 & Rlen > 2000)

svg("SB211PAT_hifiasm_canu_aln.diag.svg", width = 6, height = 6, bg=NA)
par(mfrow = c(1,1), mar = c(2.5,2.5,0,0), xpd=NA)
plot.alignments.diagonal(SB211PAT_hifiasm_canu_aln, reference_lens=SB211PAT_data$seq_len[reference_contigs], query_lens=SB211PAT_canu_data$seq_len[query_contigs],
                      reference_ori=reference_ori, query_ori=query_ori, labels_cex=0.8, angle_labels=T, labels_offset=3e5,
                      lwd=1, cols=c("#000000","#FF0000"), colour_by ="orientation", xmax=18e6, ymax=18e6)
dev.off()

################################################################################################
#####################################  plot regions only  ######################################
################################################################################################


svg("Dchry2.2_MB18102MAT_SB211PAT_chr15_synblocks_regions.svg", width = 9, height = 2, bg=NA)

par(mfrow = c(1,1), mar = c(0,0,0,0))

plot(0, cex=0, xlab="", ylab="", xaxt="n", yaxt="n", xlim = c(0,18e6), ylim = c(0,4), bty="n")


rect(1,2.8,MB18102MAT_chrom_len["chr15"], 3.2, col="gray90", border=NA)
rect(1,1.8,SB211PAT_chrom_len["chr15"], 2.2, col="gray90", border=NA)
rect(1,0.8,Dchry2.2_chrom_len["chr15"], 1.2, col="gray90", border=NA)


for (i in 1:nrow(regions)){
    
    for (j in 1:nrow(MB18102MATsyn_regions[[i]])){
        if (MB18102MATsyn_regions[[i]]$strand[j] == "+"){
            draw.horizontal.block.arrow(min(MB18102MATsyn_regions[[i]][j,c("Qstart","Qend")]), max(MB18102MATsyn_regions[[i]][j,c("Qstart","Qend")]), 3,
                                        width=0.4, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
            }
        else {
            draw.horizontal.block.arrow(max(MB18102MATsyn_regions[[i]][j,c("Qstart","Qend")]), min(MB18102MATsyn_regions[[i]][j,c("Qstart","Qend")]), 3,
                                width=0.4, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
            }
        }
    
    for (j in 1:nrow(SB211PATsyn_regions[[i]])){
        if (SB211PATsyn_regions[[i]]$strand[j] == "+"){
            draw.horizontal.block.arrow(min(SB211PATsyn_regions[[i]][j,c("Qstart","Qend")]), max(SB211PATsyn_regions[[i]][j,c("Qstart","Qend")]), 2,
                                        width=0.4, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
            }
        else {
            draw.horizontal.block.arrow(max(SB211PATsyn_regions[[i]][j,c("Qstart","Qend")]), min(SB211PATsyn_regions[[i]][j,c("Qstart","Qend")]), 2,
                                width=0.4, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
            }
        }
    
    for (j in 1:nrow(Dchry2.2syn_regions[[i]])){
        if (Dchry2.2syn_regions[[i]]$strand[j] == "+"){
            draw.horizontal.block.arrow(min(Dchry2.2syn_regions[[i]][j,c("Qstart","Qend")]), max(Dchry2.2syn_regions[[i]][j,c("Qstart","Qend")]), 1,
                                        width=0.4, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
            }
        else {
            draw.horizontal.block.arrow(max(Dchry2.2syn_regions[[i]][j,c("Qstart","Qend")]), min(Dchry2.2syn_regions[[i]][j,c("Qstart","Qend")]), 1,
                                width=0.4, width_scaler=0, length_scaler=0.1, col = region_cols[i], border=NA)
            }
        }
    
    }


rect(1,2.8,MB18102MAT_chrom_len["chr15"], 3.2, col=NA, border="gray40")
rect(1,1.8,SB211PAT_chrom_len["chr15"], 2.2, col=NA, border="gray40")
rect(1,0.8,Dchry2.2_chrom_len["chr15"], 1.2, col=NA, border="gray40")

draw.scale.bar(17e6, 3, width=1e6, height=0.2, text="1Mb", lwd=1.5, offset=0.2, cex=0.7)

dev.off()






################################################################################################
###################################  MB18102MAT to Dplex4  #####################################
################################################################################################

### this is not for main figures - just to find inversion breakpoints

MB18102MATaln <- import.paf("Dplex4_MB18102MAT.mm2asm20.paf.gz")

reference_contigs <- "chr7"
reference_ori <- c("chr7"="+")

query_contigs <- c("tig00001608", "tig00001591", "tig00001592")
query_ori <- c("tig00001608"="-", "tig00001591"="+", "tig00001592"="+")

MB18102MATaln <- subset(MB18102MATaln, Rlen>200 & reference == reference_contigs & query %in% query_contigs)

#plot.alignments.diagonal(MB18102MATaln, reference_lens=plex_data$seq_len, query_lens=MB18102MAT_data$seq_len)

plot.alignments.multi(MB18102MATaln, reference_lens=plex4_data$seq_len, query_lens=MB18102MAT_data$seq_len, reference_ori=reference_ori, query_ori=query_ori,
                      lwd=0, cols=c("#000000","#ff0000"), colour_by ="orientation", centre=FALSE, show_labels=TRUE, plot_length=18e6,
                      sigmoid=T, show_outline=FALSE, reference_above=T)



################################################################################################
#####################################  SB211PAT to Dplex4  #####################################
################################################################################################

### this is not for main figures - just to find inversion breakpoints

SB211PATaln <- import.paf("Dplex4_SB211PAT.mm2asm20.paf.gz")

reference_contigs <- "chr7"
reference_ori <- c("chr7"="+")

query_contigs <- c('h1tg000170l', 'h1tg000359l', 'h1tg000112l', 'h1tg000044l')
query_ori <- c('h1tg000170l'="-", 'h1tg000359l'="+", 'h1tg000112l'="+", 'h1tg000044l'="+")

SB211PATaln <- subset(SB211PATaln, Rlen>500 & reference == reference_contigs & query %in% query_contigs)

#plot.alignments.diagonal(SB211PATaln, reference_lens=plex_data$seq_len, query_lens=MB18102MAT_data$seq_len)

plot.alignments.multi(SB211PATaln, reference_lens=plex4_data$seq_len, query_lens=SB211PAT_data$seq_len, reference_ori=reference_ori, query_ori=query_ori,
                      lwd=0, cols=c("#000000","#ff0000"), colour_by ="orientation", centre=FALSE, show_labels=TRUE, plot_length=18e6,
                      sigmoid=T, show_outline=FALSE, reference_above=T)



SB211PATaln_sub <- subset(SB211PATaln, query == "h1tg000044l")

SB211PATaln_sub_syn <- get.synteny.blocks(SB211PATaln_sub, min_subblock_size=200)
SB211PATaln_sub_syn <- get.synteny.blocks(SB211PATaln_sub_syn, min_subblock_size=2000)
SB211PATaln_sub_syn <- get.synteny.blocks(SB211PATaln_sub_syn, min_subblock_size=20000)

plot.alignments(SB211PATaln_sub_syn, lwd=0, cols=c("#000000","#ff0000"), colour_by ="orientation", sigmoid=T, show_outline=FALSE)


plot.alignments(SB211PATaln_sub_syn, lwd=0, cols=c("#000000","#ff0000"), colour_by ="orientation", sigmoid=T, show_outline=FALSE,
                Rfirst=3.8e6, Rlast=6e6, Qfirst=1.6e6, Qlast=4.2e6)


SB211PATaln[SB211PATaln$Rstart > 3.9e6 & SB211PATaln$Rstart < 3.92e6,]

SB211PATaln[SB211PATaln$Rend > 5.93e6 & SB211PATaln$Rend < 5.99e6,]


plot.alignments(SB211PATaln_sub_syn, lwd=0, cols=c("#000000","#ff0000"), colour_by ="orientation", sigmoid=T, show_outline=FALSE,
                Rfirst=5.9e6, Rlast=8e6, Qfirst=6e6, Qlast=10e6)


par(mfrow=c(2,1))
plot.alignments(SB211PATaln_sub, lwd=0, cols=c("#000000","#ff0000"), colour_by ="orientation", sigmoid=T, show_outline=FALSE,
                Rfirst=8.8e6, Rlast=9.3e6, Qfirst=8.8e6, Qlast=9.5e6)

plot.alignments(SB211PATaln_sub_syn, lwd=0, cols=c("#000000","#ff0000"), colour_by ="orientation", sigmoid=T, show_outline=FALSE,
                Rfirst=8.8e6, Rlast=9.3e6, Qfirst=8.8e6, Qlast=9.5e6)

