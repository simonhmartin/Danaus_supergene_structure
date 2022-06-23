simple.loess.predict <- function(x, y, span, weights = NULL, max = NULL, min = NULL){
    y.loess <- loess(y ~ x, span = span, weights = weights)
    y.predict <- predict(y.loess,x)
    if (is.null(min) == FALSE) {y.predict = ifelse(y.predict > min, y.predict, min)}
    if (is.null(max) == FALSE) {y.predict = ifelse(y.predict < max, y.predict, max)}
    y.predict
    }


interleave <- function(x1,x2){
    output <- vector(length= length(x1) + length(x2))
    output[seq(1,length(output),2)] <- x1
    output[seq(2,length(output),2)] <- x2
    output
    }

setwd("~/Research/Danaus_genome/BC_REARRANGEMENTS/CNV/")

source("~/Research/dev/asynt/asynt.R")    

sample_species <- read.table("dan17.txt", as.is=T, sep="\t", header=T)


#regions (Dplex4)
regions <- as.data.frame(rbind(c(left=3909220, right=5937733),
                               c(left=5967890, right=7513824),
                               c(left=7514438, right=8804489),
                               c(left=8805142, right=9252775)))

region_cols <- c("#006437", "#a775ee", "#95e79a", "#f7614d")

### load depth data

#gene windows

# file <- "dan17.Dplex4.chr7.BT.CDS.allshared.dpstats.w1ksites.csv"
# 
# file <- "dan17.Dplex4.chr7.BT.CDS.allshared.dpstats.w5ksites.csv"
# file <- "dan17.Dplex4.chr7.BT.CDS.allshared.dpstats.w100ks10k.csv"
# 
# file <- "dan17.Dplex4.chr7.BT.CDS.allshared.dpstats.w500sitesD100k.csv"
# file <- "dan17.Dplex4.chr7.BT.CDS.allshared.dpstats.w1ksitesD100k.csv"

prefix <- "dan17.Dplex4.chr7.BT.CDS.allshared.dpstats.geneWindows"


dpstats <- read.csv(paste0(prefix,".csv"), as.is=T)
dpstats <- subset(dpstats, scaffold == "chr7" & sites >= 100)

# dpstats <- subset(dpstats, end-start+1 <= 100000)


for (ind in sample_species[,1]){
    dpstats[,paste0(ind, "_median_norm")] <- dpstats[,paste0(ind, "_median")] / median(dpstats[,paste0(ind, "_median")], na.rm=T)
    dpstats[,paste0(ind, "_mean_norm")] <- dpstats[,paste0(ind, "_mean")] / median(dpstats[,paste0(ind, "_median")], na.rm=T)
    }




species <- c("Tirumala formosa", "Danaus plexippus", "Danaus erippus", "Danaus eresimus",
                  "Danaus melanippus", "Danaus gilippus", "Danaus petilia",
                  "Danaus chrysippus chrysippus","Danaus chrysippus orientis", "Danaus chrysippus klugii")


svg(paste0(prefix, "species_hist.svg"), width = 6.5, height = 7.5, bg=NA)

par(mfrow = c(length(species),1), mar = c(0.5,1,0.1,2), bty = "n")

top=13

cols = colorRampPalette(c("black","red", "red", "red"))(top)
cols = rainbow(top, start=0.7, end=0.0)

for (spec in species){
    
    plot_inds = subset(sample_species, Species==spec)[,1]
    
    if (length(plot_inds) > 1) stat <- apply(dpstats[,paste0(plot_inds, "_median_norm")], 1, mean)
    else stat <- dpstats[,paste0(plot_inds, "_median_norm")]
    
    stat <- round(stat)
    
    stat[stat==0] <- NA
    
    col=cols[stat]
    col="black"
    
    plot(0, cex=0, ylim = c(0,top), xlim = c(0e6, 10e6), xaxt="n", yaxt="n")

    rect(0,0, 10e6,top,col="gray90", border=NA)
    
    for (i in 1:4){
        rect(regions[i,1],0, regions[i,2],top,col=paste0(region_cols[i], "55"), border=NA)
        }
    
#     segments(0, 1:top, 10e6, 1:top, col="gray80", lty=3)
    segments(0, c(1,5,10), 10e6, c(1,5,10), col="gray70", lty=3)
    
    points(dpstats$mid, stat, col=col, type="h", lwd=1)
    points(dpstats$mid, stat, col=col, type="p", pch=19, cex=0.5)
    
    axis(4,at=c(1,5,10), line=-1)

    }

dev.off()



# png("depth_and_het/Dplex4/dan21.Dplex4.chr7.BT.CDS.allshared.dpstats.geneWindows.byspecies.bars.png" , width = 2000, height = 3000, res=300, bg=NA)

species <- c("Tirumala formosa", "Danaus plexippus", "Danaus erippus", "Danaus eresimus",
                  "Danaus melanippus", "Danaus gilippus", "Danaus petilia",
                  "Danaus chrysippus chrysippus","Danaus chrysippus orientis", "Danaus chrysippus klugii")


par(mfrow = c(length(species),1), mar = c(1,1,0,2), bty = "n")

for (spec in species){
    
    plot_inds = subset(sample_species, Species==spec)[,1]
    
    if (length(plot_inds) > 1) stat <- apply(dpstats[,paste0(plot_inds, "_median_norm")], 1, mean)
    else stat <- dpstats[,paste0(plot_inds, "_median_norm")]
    
    stat <- round(stat)
#     col = cols[cut(stat, breaks = c(0, seq(0.25,10,0.5), Inf))]
    
    plot(0, cex=0, ylim = c(0,8), xlim = c(0e6, 10e6), xaxt="n", yaxt="n")

    rect(0,0, 10e6,8,col="gray90", border=NA)
    rect(regions[3,1],0, regions[3,2],8,col=paste0(region_cols[3], "55"), border=NA)

    if (spec == species[1]){
        for (r in 3){
            draw.horizontal.block.arrow(regions[r,1], regions[r,2], 6, width=2,
                                        width_scaler=0, length_scaler=0.1, col = region_cols[r], border=NA)
            }
        }


#     segments(0, 1:8, 10e6, 1:8, col="gray90", lty=3)
    rect(dpstats$start, stat, dpstats$end, 0, col="black", border=NA)
    
    lines(interleave(dpstats$start,dpstats$end) , rep(stat,each=2))
    
    axis(4,at=c(0,1,8), line=-1)

    }
# axis(1, at =seq(0,10e6,1e6), labels= seq(0,10,1))

# dev.off()






species <- c("Danaus chrysippus klugii","Danaus chrysippus orientis","Danaus chrysippus chrysippus",
             "Danaus petilia","Danaus gilippus","Danaus melanippus",
             "Danaus eresimus","Danaus erippus","Danaus plexippus")


par(mfrow = c(length(species),1), mar = c(2,2,0,1), bty = "n")

for (spec in species){
    
    plot_inds = subset(sample_species, Species==spec)[,1]
    
    if (length(plot_inds) > 1) stat <- apply(dpstats[,paste0(plot_inds, "_median_norm")], 1, mean)
    else stat <- dpstats[,paste0(plot_inds, "_median_norm")]

    col = cols[cut(stat, breaks = c(0, seq(0.25,10,0.5), Inf))]
    
    plot(0, cex=0, ylim = c(0,10), xlim = c(0e6, 10e6), xaxt="n", yaxt="n")
    abline(h=0:10, lty=3, col="gray90")
    segments(dpstats$start, stat, dpstats$end, stat, col=col, lwd=2)
    points(dpstats$start, stat, col=col, cex=1)
    points(dpstats$end, stat, col=col, cex=1)
    axis(2,at=c(0,1,5,10), line=-1)
    }
axis(1, at =seq(0,10e6,1e6), labels= seq(0,10,1))








species <- c("Danaus plexippus", "Danaus erippus", "Danaus eresimus",
                  "Danaus melanippus", "Danaus gilippus", "Danaus petilia",
                  "Danaus chrysippus chrysippus", "Danaus chrysippus orientis", "Danaus chrysippus klugii")


sp_col = c("black","gray30","gray50",
           "gray70","gray80","gray95",
           "#bc4754","#2143d1","#ffac07")

names(sp_col) <- species


inds_by_species <- sapply(species, function(species) sample_species$ID[sample_species$Species==species])
inds <- unlist(inds_by_species)

ind_col <- unlist(lapply(species, function(species) rep(sp_col[species], length(inds_by_species[[species]]))))

par(mfrow = c(1,1), mar = c(1,2,1,1), bty = "n")

fac = 6

plot(0, cex=0, ylim = c(1,fac*length(inds)+20), xlim = c(0e6, 10e6), xaxt="n")

for(i in length(inds):1){
    
    stat <- dpstats[,paste0(inds[i], "_median_norm")]
    
    idx = which(is.na(stat) ==F)
    
    x1 = dpstats$start[idx]
    x2 = dpstats$end[idx]
    y = stat[idx]
    
    rect(x1, i*fac, x2, i*fac+y, col=ind_col[i], border=NA)
    lines(interleave(x1,x2), rep(i*fac+y, each=2))
    }

axis(1, at =seq(0,10e6,1e6), labels= seq(0,10,1))









species <- c("Tirumala formosa","Danaus plexippus", "Danaus erippus", "Danaus eresimus",
                  "Danaus melanippus", "Danaus gilippus", "Danaus petilia",
                  "Danaus chrysippus chrysippus", "Danaus chrysippus orientis", "Danaus chrysippus klugii")


sp_col = c("white","black","gray30","gray50",
           "gray70","gray80","gray95",
           "#bc4754","#2143d1","#ffac07")

names(sp_col) <- species


inds_by_species <- sapply(species, function(species) sample_species$ID[sample_species$Species==species])
inds <- unlist(inds_by_species)

ind_col <- unlist(lapply(species, function(species) rep(sp_col[species], length(inds_by_species[[species]]))))

par(mfrow = c(1,1), mar = c(1,2,1,1), bty = "n")

fac = 5

plot(0, cex=0, ylim = c(1,fac*length(inds)+20), xlim = c(0e6, 10e6), xaxt="n")

for(i in length(inds):1){
    
    stat <- dpstats[,paste0(inds[i], "_median_norm")]
    stat <- ifelse(stat < 0.5, NA, stat)
    
    idx = which(is.na(stat) ==F)
    
    x = dpstats$mid[idx]
    y = stat[idx]
    w = dpstats$sites[idx]
    
    pred = simple.loess.predict(x, y, span=0.05)
    
    polygon(c(x[1], x, tail(x,1)),  c(i*fac, i*fac+pred, i*fac), col=ind_col[i])
    }

axis(1, at =seq(0,10e6,1e6), labels= seq(0,10,1))





cols = rainbow(21, start=0.6, end=0.1)

cols = colorRampPalette(c("gray90","black"))(21)

species <- c("Tirumala formosa", "Danaus plexippus", "Danaus erippus", "Danaus eresimus",
                  "Danaus melanippus", "Danaus gilippus", "Danaus petilia",
                  "Danaus chrysippus chrysippus","Danaus chrysippus orientis", "Danaus chrysippus klugii")


inds_by_species <- sapply(species, function(species) sample_species$ID[sample_species$Species==species])
inds <- unlist(inds_by_species)


par(mfrow = c(1,1), mar = c(3,2,1,1), bty = "n")

plot(0, cex=0, ylim = c(0,length(inds)), xlim = c(0e6, 10e6), xaxt="n")


rect(0, 1:length(inds)-0.3, 10e6, 1:length(inds)+0.3, col="gray90")

for (i in 1:length(inds)){
    stat <- dpstats[,paste0(inds[i], "_mean_norm")]
    
    col = cols[cut(stat, breaks = c(0, seq(0.25,10,0.5), Inf))]
    
    rect(dpstats$start, i-0.3, dpstats$end, i+0.3, col=col, border=col, lwd=1)
    
    }

axis(1, at =seq(0,10e6,1e6), labels= seq(0,10,1))


### export region 3 read depth table

columns <- paste0(sample_species[,1], "_median_norm")

#only rows in region 3, and only those at which copy number is elevated in at least two individuals
rows <- which(dpstats$start >= regions[3,1] & dpstats$start <= regions[3,2] &
              apply(dpstats[,columns], 1, function(x) sum(round(x) >=3)) > 2)

output <- round(dpstats[rows,columns])

names(output) <- sample_species[,1]

output <- cbind(dpstats[rows, c("start","end")], output)

write.table(output, paste0(prefix, "Region3_highCopy.tsv"), quote=F, row.names=F)

