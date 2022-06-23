
get.chrom.pos <- function(scaffold, position, scaf_len, scaf_ori){
    
    npos <- length(position)
    
    offsets <- cumsum(scaf_len) - scaf_len
    
    chrom_pos <- numeric(length = npos)
    
    for (i in 1:npos){
        if (scaf_ori[scaffold[i]] == "+") chrom_pos[i] <- offsets[scaffold[i]] + position[i]
        else chrom_pos[i] <- offsets[scaffold[i]] + scaf_len[scaffold[i]] - position[i]
        }
    
    chrom_pos
    }

get.contig.lengths <- function(fai_file){
    fai <- read.table(fai_file, as.is=T)
    lengths <- fai[,2]
    names(lengths) <- as.character(fai[,1])
    lengths
    }

###################################################################################################

dir="/home/simon/Research/Danaus_genome/BC_REARRANGEMENTS/depth/"

setwd(dir)



# genomes <- c("Dchry2.2", "Dchry2HAP", "MB18102PAT", "MB18102MAT", "SB211PAT", "SB211MAT")

genomes <- c("Dchry2.2", "MB18102MAT", "SB211PAT")


fai_files <- c(Dchry2.2="../../Dchry2/Dchry2.2.fa.fai",
               Dchry2HAP="../../Dchry2/Dchry2.haplotigs.fasta.fai",
               MB18102MAT="../../MB181_trio/MB18102_MAT.fasta.fai",
               MB18102PAT="../../MB181_trio/MB18102_PAT.fasta.fai",
               SB211MAT="../../SB211_trio/SB211_MAT.fasta.fai",
               SB211PAT="../../SB211_trio/SB211_PAT.fasta.fai")

contig_len <- sapply(fai_files, get.contig.lengths, simplify=FALSE)

contig_ori <- list(Dchry2.2=c("contig15.1"="+"),
                   Dchry2HAP=c("000397F"="-"),
                   SB211PAT=c('h1tg000170l'="-",  'h1tg000359l'="+", 'h1tg000112l'="+", 'h1tg000044l'="+"),
                   SB211MAT=c('tig00000001'="+",  'tig00000185'="-"),
                   MB18102PAT=c('tig00000001'="-"),
                   MB18102MAT=c("tig00001608"="-", "tig00001591"="+", "tig00001592"="+"))

target_contigs <- sapply(contig_ori, names, simplify=F)

contig_len <- sapply(genomes, function(genome) contig_len[[genome]][target_contigs[[genome]]])

### divergence data
dp_data <- sapply(genomes, function(genome) read.csv(paste0("chry10.BT.", genome, ".DPstats.w50.csv"), as.is = T), simplify=FALSE)


### population data names of the labelling populations

samples <- list(Dchry2.2=c("SM15W61", "SM15W66", "SM15W72", "SM15W74"),
                MB18102PAT=c("SM15W61", "SM15W66", "SM15W72", "SM15W74"),
                MB18102MAT=c("SM16N04", "SM16N05", "SM16N20", "SM16N37"),
                Dchry2HAP=c("SM16N04", "SM16N05", "SM16N20", "SM16N37"),
                SB211PAT=c("SM16S06", "SM17S01")
                )

### add chromosome data

for (genome in genomes){
    #adjust positions for chromosomal plot
    for (header in c("start","end","mid")) {
        dp_data[[genome]][,header] <- get.chrom.pos(dp_data[[genome]]$scaffold, dp_data[[genome]][,header],
                                                     contig_len[[genome]], contig_ori[[genome]])
        }
    }


### compute median

for (genome in genomes){
    for (s in samples[[genome]]) {
        dp_data[[genome]][,paste0(s, "_median_norm")] <- dp_data[[genome]][,paste0(s, "_median")] / median(dp_data[[genome]][,paste0(s, "_median")], na.rm=T)
        }
    
    dp_data[[genome]]$median_norm <- apply(dp_data[[genome]][,paste0(samples[[genome]], "_median_norm")], 1, mean, na.rm=T)
    
    }


### plot dp with best match

labels <- c(Dchry2.2="D. chrysippus klugii allele (Dchry2.2 assembly)",
            MB18102PAT="D. chrysippus klugii allele (MB18102PAT assembly)",
            MB18102MAT="D. chrysippus chrysippus allele (MB18102MAT assembly)",
            Dchry2HAP="D. chrysippus chrysippus allele (Dchry2 alternative haplotig)",
            SB211PAT="D. chrysippus orientis allele (SB211PAT assembly)")


svg(paste0(paste(genomes, collapse="_"), ".DP.w50.svg"), width=12, height=8, bg=NA)

par(mfrow=c(length(genomes), 1), mar = c(4,4,2,1), bty = "n")


for (genome in genomes){
    
    plot(dp_data[[genome]]$mid, dp_data[[genome]]$median_norm, xlim = c(0,18e6), ylim = c(0,3), xlab="", ylab = "", xaxt="n", yaxt="n")
    
    segments(0,1,sum(contig_len[[genome]]), 1, col="blue", lty=3)
    
    rect(cumsum(contig_len[[genome]]), 0, cumsum(contig_len[[genome]])-contig_len[[genome]], 3)
    
    xmax = ceiling(sum(contig_len[[genome]])/1e6)
    
    axis(1, at = seq(0,xmax*1e6,1e6), labels=0:xmax)
    axis(2, line=-3)
    
    mtext(side=1,"Position on chromosome 15 (Mb)", line=3)
    mtext(side=2, "Normalised read depth", line=0)
    
    mtext(side=3, at=2e5, text=labels[[genome]], adj=0, line=0)
    
    }

dev.off()

