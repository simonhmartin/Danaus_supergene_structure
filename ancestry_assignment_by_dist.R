
which.min.threshold <- function(values, threshold = 1, mode="ratio", nomin=0){
    sorted_values <- sort(values)
    #if all but one is NA, then return the minimum
    if (length(sorted_values) == 1) return(which.min(values))
    
    if (mode == "ratio"){
        try(if (sorted_values[1]/sorted_values[2] <= threshold) return(which.min(values))
            else return(nomin), silent=TRUE)
        }
    else {
        try(if (sorted_values[2] - sorted_values[1] >= threshold) return(which.min(values))
            else return(nomin), silent=TRUE)
        }
    NA
    }

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

dir="/home/simon/Research/Danaus_genome/BC_REARRANGEMENTS/short_read_analyses/"

setwd(dir)


###########################################################################################################
##########################################  individual plots  #############################################
###########################################################################################################


genomes <- c("Dchry2.2", "Dchry2HAP", "MB18102PAT", "MB18102MAT", "SB211PAT", "SB211MAT")

prefixes <- paste0("chry10.BT.", genomes, ".DP8GQ20.chr15.divStats.w25ksitesMin10k")
names(prefixes) <- genomes

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
div_data <- sapply(genomes, function(genome) read.csv(paste0(dir,prefixes[genome],".csv.gz"), as.is = T), simplify=FALSE)


### population data names of the labelling populations
pop_names <- c("klugii", "chrysippus", "orientis")

pop_cols <- c("#ffac07","#bc4754","#2143d1","gray90")
names(pop_cols) <- pop_names


### add chromosome data

for (genome in genomes){
    #remove very long windows
    too_long <- div_data[[genome]]$end - div_data[[genome]]$start +1 > 1e5
    div_data[[genome]][too_long,-1:-5] <- NA
    #adjust positions for chromosomal plot
    for (header in c("start","end","mid")) {
        div_data[[genome]][,header] <- get.chrom.pos(div_data[[genome]]$scaffold, div_data[[genome]][,header],
                                                     contig_len[[genome]], contig_ori[[genome]])
        }
    }


### compute best match

for (genome in genomes){
    #da divergence between reference and each population
    for (pop in pop_names){
        div_data[[genome]][,paste("da", pop, "REF", sep="_")] <- div_data[[genome]][,paste("dxy", pop, "REF", sep="_")] - div_data[[genome]][,paste0("pi_", pop)]
        }
    
    div_data[[genome]] <- div_data[[genome]][order(div_data[[genome]]$start),]
    
    div_data[[genome]]$best_match <- apply(div_data[[genome]][,paste("da", pop_names, "REF", sep="_")], 1, which.min.threshold,
                                           threshold=0.01, mode="difference", nomin=4)
    
    }


### plot divergence with best match

svg("all_chr15.best_match_da0.01.overlay.svg", width=12, height=8, bg=NA)

layout(rbind(matrix(1:4, ncol=2), matrix(5:8, ncol=2), matrix(9:12, ncol=2)),
       heights=c(3,1,3,1,3,1))

for (genome in genomes){
        
#     layout(matrix(c(1,2), ncol=1), heights=c(3,1))
    par(mar = c(0,4,2,0), bty = "n")
    
    box_left = 0
    box_right = max(div_data[[genome]]$end)/1e6
    xlim = c(0,18)
    plot(0, cex=0, xlim = xlim, ylim = c(-0.01,0.08), ylab = "", xaxt="n")
    for (pop in pop_names){
#         lines(div_data[[genome]]$mid/1e6, div_data[[genome]][,paste("da", pop, "REF", sep="_")],
#               type = "l", col = pop_cols[pop], lwd=2)
        points(div_data[[genome]]$mid/1e6, div_data[[genome]][,paste("da", pop, "REF", sep="_")],
              type = "p", col = pop_cols[pop], pch=19, cex=.5)
        }
    mtext(2, text=expression(italic("d"["a"])), line=2, cex=1)
    mtext(3, at=0, text=genome, adj=0, line=-2)
    
    par(mar = c(2,4,0,0))
    width <- 0.8
    plot(0, cex = 0, xlim = xlim, ylim = c(0,1.6), ylab = "", xaxt="n", yaxt = "n", xlab = "")
    rect(div_data[[genome]]$start/1e6, 1+width/2, div_data[[genome]]$end/1e6, 1-width/2,
         col = pop_cols[div_data[[genome]]$best_match], border=NA, lwd=0)
    rect(box_left, 1-width/2, box_right, 1+width/2)
    axis(1, at=0:18, labels=0:18)
    }

dev.off()


### plot best match only for inclusion in alignment figures

for (genome in genomes){
    svg(paste0(prefixes[[genome]],".best_match_da0.01.overlay.svg"), width = 9, height = 4.5, bg=NA)
    par(mfrow = c(1,1), mar = c(0,2.5,0,0))
    plot(0,cex = 0, xlim = c(1, 18e6), ylim = c(1.3,-.3), xlab = "", ylab = "", bty = "n", yaxt="n", xaxt="n")
    rect(div_data[[genome]]$start, 1, div_data[[genome]]$end, 1.1, col = pop_cols[div_data[[genome]]$best_match], border=NA, lwd=0)
    box_right = sum(contig_len[[genome]])
    rect(0, 1, box_right, 1.1, border="gray40")
    dev.off()
    }

### and for pairs

for (pair in list(c("MB18102PAT","MB18102MAT"), c("Dchry2.2","MB18102MAT"), c("MB18102MAT","SB211PAT"), c("SB211PAT","MB18102PAT"), c("SB211PAT","Dchry2.2"))){
    svg(paste0(prefixes[[pair[1]]], "_", prefixes[[pair[2]]], ".best_match_da0.01.overlay.svg"), width = 7, height = 3.5, bg=NA)
    par(mfrow = c(1,1), mar = c(0,2.5,0,0))
    plot(0,cex = 0, xlim = c(1, 18e6), ylim = c(1.3,-.3), xlab = "", ylab = "", bty = "n", yaxt="n", xaxt="n")
    
    rect(div_data[[pair[1]]]$start, 0, div_data[[pair[1]]]$end, -.1, col = pop_cols[div_data[[pair[1]]]$best_match], border=NA, lwd=0)
    rect(0, 0, sum(contig_len[[pair[1]]]), -.1, border="gray40")

    rect(div_data[[pair[2]]]$start, 1, div_data[[pair[2]]]$end, 1.1, col = pop_cols[div_data[[pair[2]]]$best_match], border=NA, lwd=0)
    rect(0, 1, sum(contig_len[[pair[2]]]), 1.1, border="gray40")
    dev.off()
    }



