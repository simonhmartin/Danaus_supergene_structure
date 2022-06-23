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


genomes <- c("Dchry2.2", "MB18102PAT", "MB18102MAT", "SB211PAT", "SB211MAT")


fai_files <- c(Dchry2.2="../../Dchry2/Dchry2.2.fa.fai",
               MB18102MAT="../../MB181_trio/MB18102_MAT.fasta.fai",
               MB18102PAT="../../MB181_trio/MB18102_PAT.fasta.fai",
               SB211MAT="../../SB211_trio/SB211_MAT.fasta.fai",
               SB211PAT="../../SB211_trio/SB211_PAT.fasta.fai")

contig_len <- sapply(fai_files, get.contig.lengths, simplify=FALSE)

contig_ori <- list(Dchry2.2=c("contig15.1"="+"),
                   SB211PAT=c('h1tg000170l'="-",  'h1tg000359l'="+", 'h1tg000112l'="+", 'h1tg000044l'="+"),
                   SB211MAT=c('tig00000001'="+",  'tig00000185'="-"),
                   MB18102PAT=c('tig00000001'="-"),
                   MB18102MAT=c("tig00001608"="-", "tig00001591"="+", "tig00001592"="+"))

target_contigs <- sapply(contig_ori, names, simplify=F)

contig_len <- sapply(genomes, function(genome) contig_len[[genome]][target_contigs[[genome]]])

#regions
regions <- sapply(genomes, function(g) read.table(paste0("../alignment_plots/", g,"_region_coordinates.tsv"),
                                                  header=T, as.is=T), simplify=F)

region_cols <- c("#006437", "#a775ee", "#95e79a", "#f7614d")


### population data names of the labelling populations
pop_names <- c("klugii", "chrysippus", "orientis")

pop_cols <- c("#ffac07","#bc4754","#2143d1","gray90")
names(pop_cols) <- pop_names


### divergence data
prefixes <- paste0("chry10.BT.", genomes, ".DP8GQ20.chr15.divStats.w25ksitesMin20k")
names(prefixes) <- genomes

div_data <- sapply(genomes, function(genome) read.csv(paste0(dir,prefixes[genome],".csv.gz"), as.is = T), simplify=FALSE)


### add chrom pos and da for pairs
for (genome in genomes){
    #remove very long windows
    too_long <- div_data[[genome]]$end - div_data[[genome]]$start +1 > 1e5
    div_data[[genome]][too_long,-1:-5] <- NA

    #adjust positions for chromosomal plot
    for (header in c("start","end","mid")) {
        div_data[[genome]][,header] <- get.chrom.pos(div_data[[genome]]$scaffold, div_data[[genome]][,header],
                                                     contig_len[[genome]], contig_ori[[genome]])
        }
    
    #da divergence between reference and each population
    for (i in 1:(length(pop_names)-1)){
        pop1 <- pop_names[i]
        for (j in (i+1):length(pop_names)){
            pop2 <- pop_names[j]
            div_data[[genome]][,paste("da", pop1, pop2, sep="_")] <- (div_data[[genome]][,paste("dxy", pop1, pop2, sep="_")] -
                                                                      0.5*(div_data[[genome]][,paste0("pi_", pop1)] + div_data[[genome]][,paste0("pi_", pop2)]))
            }
        }
    
    div_data[[genome]] <- div_data[[genome]][order(div_data[[genome]]$start),]
    
    }

### add mean stats for regions
for (genome in genomes){
    for (i in 1:nrow(regions[[genome]])){
        windows <- which(div_data[[genome]]$start >= regions[[genome]][i, "Qstart"] &
                         div_data[[genome]]$end <= regions[[genome]][i, "Qend"])
        
        for (stat in names(div_data[[1]])[-(1:5)]){
            regions[[genome]][i,stat] <- mean(div_data[[genome]][windows,stat], na.rm=T)
            }
        }
    }
            


comps <- list(c("klugii", "chrysippus"),
              c("chrysippus", "orientis"),
              c("klugii", "orientis"))


stat = "da"
var_label <- expression(italic("d"["a"]))
ylim <- c(-0.005,0.068)

# stat = "dxy"
# var_label <- expression(italic("d"["XY"]))
# ylim <- c(0,0.09)



xlim = c(1, 18e6)


for (reference in genomes){
    for (comp in comps) {
        varname = paste(stat, comp[1], comp[2], sep="_")
            
            svg(paste0(prefixes[reference],".",comp[1], "_", comp[2], ".", stat, ".svg"), width = 7, height = 3, bg=NA)
            par(mar = c(3,2.5,0,0), bty = "n")
            plot(div_data[[reference]]$mid, div_data[[reference]][,varname], xlim = xlim, ylim = ylim,
                xaxt="n", yaxt="n", xlab = "", ylab = "", bty = "n", pch=19, cex=0.7, type="p")
            axis(1, at=seq(0, xlim[2], 1e6), labels = 0:(xlim[2]/1e6))
            axis(2, line=-1.2, )
            mtext(2, text=var_label, line=1, cex=1.5)

            #add region mean
            # for (i in 1:nrow(regions[[reference]])){
            #     segments(regions[[reference]][i,"Qstart"], regions[[reference]][i,varname],
            #              regions[[reference]][i,"Qend"], regions[[reference]][i,varname],
            #              col=region_cols[regions[[reference]][i,"Region"]])
            #     }

            rect(regions[[reference]][,"Qstart"], 0, regions[[reference]][,"Qend"], ylim[2],
                 col=paste0(region_cols[regions[[reference]][,"Region"]], "33"), border=NA)

            dev.off()
        }
    }

