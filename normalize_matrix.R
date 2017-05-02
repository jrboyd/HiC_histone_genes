library(Sushi)
library(GenomicRanges)
library(magrittr)
source("functions_my_HiC_plots.R")
files = dir("bins_balanced_filtered2//", full.names = T)
# files = files[2]
# files = files[c(1,3:4)]
if(file.exists("hist_ref.save")){
  load("hist_ref.save")
}else{
  ref = parse_gtf("gencode.v24.annotation.gtf", rownames_attrib = "gene_id", feature_type = "gene")
  is_hist = grepl("^HIST", ref$gene_name)
  hist_ref = ref[is_hist,]
  ref[ref$gene_name == "NPAT",]
  hist_ref = rbind(hist_ref, ref[ref$gene_name == "NPAT",], ref[ref$gene_name == "HINFP" ,])
  save(hist_ref, file = "hist_ref.save")
}
loci_gap = 10^6
hist_pos = GRanges(seqnames = hist_ref$chrm, IRanges(hist_ref$start - loci_gap, hist_ref$end + loci_gap))
names(hist_pos) = hist_ref$gene_id
hist_loci = reduce(hist_pos)

all_chrm = unique(sapply(strsplit(rownames(binned), " "), function(x)x[1]))
names(all_chrm) = all_chrm
chrm_sizes = sapply(all_chrm, function(x)grepl(paste0(x, " "), rownames(binned)) %>% sum)

all_binned = list()
pdf("matrix_normalizing2.pdf")
for(f in files){
  key = basename(f) %>% sub(pattern = ".save", replacement = "")
  print(paste(key, ":", f))
  outf = paste0("bins_balanced_filtered2_normalized/", key, ".save")
  if(file.exists(outf)){
    print("skip, outf exists")
    next
  }
  load(f)
  main = mean(diag(binned))
  print("byChr")
  binned_byChr = lapply(all_chrm, function(x){
    k = grepl(paste0(x, " "), rownames(binned))
    binned[k,k]
  })
  print("antiChr")
  binned_antiChr = lapply(all_chrm, function(x){
    k = grepl(paste0(x, " "), rownames(binned))
    binned[k,!k]
  })
  
  others = sapply(1:20, function(x){
    print(x)
    sapply(binned_byChr, function(y){
      mean(diag(y[,-1:-x]))
    })
  })
  distant = sapply(801:820, function(x){
    print(x)
    sapply(binned_byChr, function(y){
      mean(diag(y[,-1:-x]))
    })
  })
  distant_means = apply(distant, 1, mean)
  
  inter = sapply(binned_antiChr, function(x){
    set.seed(0)
    mean(x)
    
  })
  plot(distant_means, inter, type = "n")
  text(distant_means, inter, labels = names(distant_means))
  title(key)
  
  metrics = list(main = main, near = others, distant_means = distant_means, inter = inter)
  save(metrics, bin_size, file = outf)
  # all_binned[[key]] = binned
  
}
dev.off()

files = dir("bins_balanced_filtered2_normalized/", full.names = T)
m = matrix(0, ncol = 0, nrow = 23)
for(f in files){
  load(f)
  print(f)
  m = cbind(m ,metrics$distant_means / metrics$inter)
}
colnames(m) = basename(files)
library(gplots)
m = apply(m, 2, function(x)x/max(x))
heatmap.2(m, Colv = T, Rowv = F, trace = "n", margins = c(10,5), cexCol = 1)

x = distant
plot(0, xlim = c(1,20), ylim = range(distant))
distant_means = apply(distant, 1, mean)
boxplot(distant_means / min(distant_means))
# 
# loci_gap = 10^6
# hist_pos = GRanges(seqnames = hist_ref$chrm, IRanges(hist_ref$start - loci_gap, hist_ref$end + loci_gap))
# names(hist_pos) = hist_ref$gene_id
# hist_loci = reduce(hist_pos)
# 
# #p_ext - how far plot is extended around loci
# my_plotHic = function(binned_name, loci_gr, gene_ref, p_ext = 5*10^5){
#   binned = all_binned[[binned_name]]
#   
#   # names(loci) = c("HIST2.1", "HIST2.2", "HIST2.3")
#   for(i in 1:length(loci_gr)){
#     #range of interest
#     roi = list(chrm = as.character(seqnames(loci_gr))[i], start = start(loci_gr)[i], end = end(loci_gr)[i])
#     roi.plotted = roi #range to be plotted
#     
#     roi.plotted$start = max(roi.plotted$start - p_ext, 0)
#     roi.plotted$end = roi.plotted$end + p_ext
#     # roi.plotted$start = 25000000
#     # roi.plotted$end = 29000000
#     
#     k = grepl(roi.plotted$chrm, rownames(binned))
#     m = binned[k,k]
#     s = as.integer(roi.plotted$start / bin_size)
#     e = as.integer(roi.plotted$end / bin_size) + 1
#     m = m[s:e, s:e]
#     rownames(m) = (s:e - 1) * bin_size
#     colnames(m) = rownames(m)
#     
#     roi.plotted.gr = GRanges(seqnames = roi.plotted$chrm, IRanges(roi.plotted$start, roi.plotted$end))
#     olaps = findOverlaps(query = roi.plotted.gr, subject = hist_pos)
#     roi.genes = names(hist_pos[subjectHits(olaps)])
#     roi.genes.sym = gene_ref[roi.genes,]$gene_name
#     roi.genes.pos = rowMeans(gene_ref[roi.genes,4:5])
#     
#     colors = function(n)rgb(colorRamp(c("white", "blue", "red"))((1:n-1)/n)/255)
#     phic = plotHic(log10(as.matrix(m*10000)+1), chrom = roi.plotted$chrm, chromstart = roi.plotted$start, chromend = roi.plotted$end, palette = colors, zrange = c(0,2.5), max_y = sqrt(nrow(m)^2/2))
#     title(binned_name)
#     axis(side = 1, at = roi.genes.pos, labels = roi.genes.sym, line = 1)
#     labelgenome(roi$chrm,roi.plotted$start,roi.plotted$end,side=1,scipen=20,n=4,scale="Mb",edgeblankfraction=0.20,chromline=.5,scaleline=0.5, line = 3)
#     addlegend(phic[[1]],palette=phic[[2]],title="score",side="right",bottominset=0.4,topinset=0,xoffset=-.035,labelside="left",width=0.025,title.offset=0.035)
#   }
# }
# 
# setwd("output_hmaps")
# ext = 3*10^6
# for(i in 1:length(all_binned)){
#   for(j in 1:length(hist_loci)){
#     root = paste(sep = "_",
#                  as.character(seqnames(hist_loci)[j]),
#                  paste0(as.integer(mean(c(start(hist_loci)[j], end(hist_loci)[j]))/10^3), "K"),
#                  names(all_binned)[i])
#     print(root)
#     # png(paste0(root, ".png"), width = 800, height = 600)
#     # my_plotHic(names(all_binned)[i], loci_gr = hist_loci[j], hist_ref, p_ext = ext)
#     # # chrm = as.character(seqnames(hist_loci)[j])
#     # # start = start(hist_loci)[j]
#     # # end = end(hist_loci)[j]
#     # # my_plotHic_lines(binned_name = names(all_binned)[i], 
#     # #                  chrm = chrm, start = start - ext, end = end + ext, 
#     # #                  gene_ref = hist_ref, 
#     # #                  min_dist = 4, 
#     # #                  n_displayed = 200)
#     # dev.off()
#   }
#   loci = hist_loci[11]
#   end(loci) = end(hist_loci[12])
#   png(paste0("NPAT-and-HINFP_", names(all_binned)[i], ".png"))
#   my_plotHic(binned_name = names(all_binned)[i], loci_gr = loci, gene_ref = hist_ref, p_ext = ext)
#   dev.off()
#   # chrm = as.character(seqnames(loci))
#   # start = start(loci)
#   # end = end(loci)
#   # my_plotHic_lines(binned_name = names(all_binned)[i],
#   #                  chrm = chrm, start = start - ext, end = end + ext,
#   #                  gene_ref = hist_ref,
#   #                  min_dist = 10,
#   #                  n_displayed = 200)
#   # dev.off()
# }
