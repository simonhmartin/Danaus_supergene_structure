source("~/Research/dev/asynt/asynt.R")

setwd("/home/simon/Research/Danaus_genome/BC_REARRANGEMENTS/TE_content/")

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


############################################

genomes=c("DplexMex","Dchry2.2", "MB18102MAT","SB211PAT")


fai_files <- c(DplexMex="../../Dplex_mex/dplex_mex.fa.fai",
               Dchry2.2="../../Dchry2/Dchry2.2.fa.fai",
               MB18102MAT="../../MB181_trio/MB18102_MAT.fasta.fai",
               SB211PAT="../../SB211_trio/SB211_PAT.fasta.fai")

contig_len <- sapply(fai_files, get.contig.lengths, simplify=FALSE)

contig_ori <- list(DplexMex=c("mxdp_6"="-"),
                   Dchry2.2=c("contig15.1"="+"),
                   SB211PAT=c('h1tg000170l'="-",  'h1tg000359l'="+", 'h1tg000112l'="+", 'h1tg000044l'="+"),
                   MB18102MAT=c("tig00001608"="-", "tig00001591"="+", "tig00001592"="+"))

target_contigs <- sapply(contig_ori, names, simplify=F)

contig_len <- sapply(genomes, function(genome) contig_len[[genome]][target_contigs[[genome]]])


############ TE data ###################

types <- c("DNA",
"Helitron",
"LINE",
"LTR",
"Retroposon",
"rRNA",
"Satellite",
"SINE",
"tRNA",
"Unknown")

te_data <- sapply(genomes, function(g) read.table(paste0(g, ".chr15.TE_50kb"), sep= "\t", as.is=T), simplify=F)


#add start and convert to chrom pos
for (genome in genomes){
    names(te_data[[genome]]) <- c("scaf","end",types)
    
    te_data[[genome]]$start <- te_data[[genome]]$end - 49999
    
    #adjust positions for chromosomal plot
    for (header in c("start","end")) {
        te_data[[genome]][,header] <- get.chrom.pos(te_data[[genome]]$scaf, te_data[[genome]][,header],
                                                     contig_len[[genome]], contig_ori[[genome]])
        }
    
    te_data[[genome]] <- te_data[[genome]][order(te_data[[genome]]$start),]
    
    for (type in types) te_data[[genome]][is.na(te_data[[genome]][,type])==TRUE,type] <- 0
    }


#set lower and upper bounds for rectangles
te_upper <- list()
te_lower <- list()

for (genome in genomes){
    te_upper[[genome]]$start <- te_data[[genome]]$end - 49999
    
    te_upper[[genome]] <- t(apply(te_data[[genome]][,types], 1, function(x) cumsum(x)))
    te_lower[[genome]] <- te_upper[[genome]] - te_data[[genome]][,types]
    }
    

cols <- c(
"#0075DC", #Blue
"#00998F", #Turquoise
"#FFA405", #Orpiment
"#FF0010", #Red
"#740AFF", #Violet
"#426600", #Quagmire
"#C20088", #Mallow
"#9DCC00", #Lime
"#F0A3FF", #Amethyst
"gray80")

# "#003380", #Navy
# "#94FFB5") #Jade

# cols = c(rainbow(9), "gray80")
names(cols) <- types

############# region data ################

regions <- sapply(genomes, function(g) read.table(paste0("../alignment_plots/", g,"_region_coordinates.tsv"),
                                                  header=T, as.is=T), simplify=F)

for (genome in genomes) {
    names(regions[[genome]])[names(regions[[genome]]) == "Qstart"] <- "start"
    names(regions[[genome]])[names(regions[[genome]]) == "Qend"] <- "end"
    }

region_cols <- c("#006437", "#a775ee", "#95e79a", "#f7614d")

####################################################


labels <- c(DplexMex="D. plexippus (MEX_DaPlex assembly)",
               Dchry2.2="D. chrysippus klugii allele (Dchry2.2 assembly)",
               MB18102MAT="D. chrysippus chrysippus allele (MB18102MAT assembly)",
               SB211PAT="D. chrysippus orientis allele (SB211PAT assembly)")




svg("DplexMex_Dchry2.2_MB18102MAT_SB211PAT.TE_50kb.svg", width=12, height=12)

par(mar=c(4,2,4,1), mfrow=c(length(genomes), 1))


for (genome in genomes){
    plot(0, cex=0, xlim = c(0, 18e6), ylim = c(-5000, 50000), bty = "n", xaxt="n", yaxt="n", xlab="", ylab="")

    for (type in types){
        rect(te_data[[genome]]$start, te_lower[[genome]][,type], te_data[[genome]]$end, te_upper[[genome]][,type], border=NA, col=cols[type])
        }


    for (i in 1:nrow(regions[[genome]])){
        if (regions[[genome]]$strand[i] == "+"){
            draw.horizontal.block.arrow(min(regions[[genome]][i,c("start","end")]), max(regions[[genome]][i,c("start","end")]), -3000,
                                            width=4000, width_scaler=0, length_scaler=0.1, col = region_cols[regions[[genome]][i,"Region"]], border=NA)
            }
        else {
            draw.horizontal.block.arrow(max(regions[[genome]][i,c("start","end")]), min(regions[[genome]][i,c("start","end")]), -3000,
                                width=4000, width_scaler=0, length_scaler=0.1, col = region_cols[regions[[genome]][i,"Region"]], border=NA)
            }
        }
    
    xmax = ceiling(sum(contig_len[[genome]])/1e6)
    
    axis(1, at = seq(0,xmax*1e6,1e6), labels=0:xmax)
    axis(2, line=-3)

    mtext(side=1,"Position on chromosome 15 (Mb)", line=3)
    mtext(side=2, "Repeat content (bp)", line=0)
    
    mtext(side=3, at=2e5, text=labels[[genome]], adj=0, line=-1)
    
    if (genome == genomes[1]) legend("topright", fill=cols, legend=types, border=NA, ncol=2)

    }

dev.off()


### wider

svg("DplexMex_Dchry2.2_MB18102MAT_SB211PAT.TE_50kb.wide.svg", width=15, height=10, bg=NA)

par(mar=c(4,2,1,0), mfrow=c(length(genomes), 1))


for (genome in genomes){
    plot(0, cex=0, xlim = c(0, 18e6), ylim = c(-5000, 50000), bty = "n", xaxt="n", yaxt="n", xlab="", ylab="")

    for (type in types){
        rect(te_data[[genome]]$start, te_lower[[genome]][,type], te_data[[genome]]$end, te_upper[[genome]][,type], border=NA, col=cols[type])
        }


    for (i in 1:nrow(regions[[genome]])){
        if (regions[[genome]]$strand[i] == "+"){
            draw.horizontal.block.arrow(min(regions[[genome]][i,c("start","end")]), max(regions[[genome]][i,c("start","end")]), -3000,
                                            width=4000, width_scaler=0, length_scaler=0.1, col = region_cols[regions[[genome]][i,"Region"]], border=NA)
            }
        else {
            draw.horizontal.block.arrow(max(regions[[genome]][i,c("start","end")]), min(regions[[genome]][i,c("start","end")]), -3000,
                                width=4000, width_scaler=0, length_scaler=0.1, col = region_cols[regions[[genome]][i,"Region"]], border=NA)
            }
        }
    
    xmax = ceiling(sum(contig_len[[genome]])/1e6)
    
    axis(1, at = seq(0,xmax*1e6,1e6), labels=0:xmax)
    axis(2, line=-3)

    if (genome == genomes[4]) mtext(side=1,"Position on chromosome 15 (Mb)", line=3)
    mtext(side=2, "Repeat content (bp)", line=0)
    
    #mtext(side=3, at=2e5, text=labels[[genome]], adj=0, line=-1)
    
    if (genome == genomes[1]) legend("topright", fill=cols, legend=types, border=NA, ncol=3)

    }

dev.off()














