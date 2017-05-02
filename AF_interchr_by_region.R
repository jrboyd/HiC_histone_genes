#summarizes interchromosomal interactions focused on a particular region
library(magrittr)
library(GenomicRanges)
source("functions_my_HiC_plots.R")

ref = parse_gtf("gencode.v24.annotation.gtf", rownames_attrib = "gene_id", feature_type = "gene")
goi = c("ZEB1", "PLAUR", "IL6", "CDH1", "FN1")
# goi = "FN1"
goi_ref = ref[0,]
for(g in goi){
  goi_ref = rbind(goi_ref,subset(ref, gene_name == g))
}

goi_loci = GRanges(seqnames = goi_ref$chrm, IRanges(goi_ref$start, goi_ref$end), gene_name = goi_ref$gene_name)




files = dir("bins_balanced_filtered3/", full.names = T, pattern = "MCF10")
files= files[c(3,1,2)]

all_binned = list()
for(f in files){
  key = basename(f) %>% sub(pattern = ".save", replacement = "")
  print(f)
  load(f)
  all_binned[[key]] = binned
}

plot_region_interchromosomal = function(binned, bin_size, main, chrm, start, end, mask_n, best_n){
  extract = extract_HiC_chrmRange(M = binned, chrm = chrm, start = start, end = end, bin_size = bin_size, by = "r")
  extract_means =  apply(extract, 2, mean)
  to_mask = genomic2matrix_index(start- bin_size * mask_n, end + bin_size * mask_n, bin_size)
  kept = setdiff(names(extract_means), paste(chrm, to_mask))
  extract_means = extract_means[kept]
  o = order(extract_means, decreasing = T)
  
  # best =  head(extract_means[o], n = best_n)
  best = extract_means[extract_means>0]
  norm_best = best / max(best)
  # plot(norm_best[1:1000])
  # for(i in 1:(best_n-10)){
  #   k = i:(i+10)
  #   range(best[k])
  # }
  n = length(best)
  mid = mean(c(start, end))
  i = genomic2matrix_index(mid, mid, bin_size)
  y1 = which(rownames(binned) == paste(chrm, i))
  x1s = rep(0, n)
  y1s = rep(y1, n)
  y2s = sapply(names(best), function(x)which(rownames(binned)==x))
  x2s = rep(1, n)
  
  par(mai = rep(1,4))
  plot(0, xlim = c(-.05, 1.05), ylim = c(1, nrow(binned)), type = "n", axes = F, ylab = "", xlab = "")
  
  for(i in n:1){
    weight = 1-norm_best[i]
    y = y2s[i]
    lines(c(0,1), c(y1, y), lwd  = 2, col = rgb(weight,weight,weight))
  }
  
  
  
  chrm_starts = cumsum(c(0, chrm_sizes[-length(chrm_sizes)]))+1
  chrm_ends = cumsum(chrm_sizes)
  chrm_mids = rowMeans(cbind(chrm_starts, chrm_ends))
  for(i in 1:length(chrm_sizes)){
    s = chrm_starts[i]
    e = chrm_ends[i]
    rect(xleft = -.05, ybottom = e, xright = -.02, ytop = s)
    rect(xleft = 1.02, ybottom = e, xright = 1.05, ytop = s)
  }
  title(main = main, sub = "Intrachrosomal interactions ignored.\nTop 1000 interactions shown.")
  legend("bottom", legend = "highest relative frequency", col = "black", lty = 1, lwd = 2)
  axis(side = 2, at = chrm_mids, labels = names(chrm_sizes), las = 2, tick = F, cex.axis = .6)
  axis(side = 4, at = chrm_mids, labels = names(chrm_sizes), las = 2, tick = F, cex.axis = .6)
  names(best) = sapply(strsplit(names(best), " "), function(x){
    paste(x[1], (as.integer(x[2]) - 1) * bin_size)
  })
  return(best)
}

outdir = "output_AF_hg38_integrated_filtered3_interchromosomal_v2/"
dir.create(outdir, showWarnings = F)
outdir = normalizePath(outdir)
setwd(outdir)

for(j in (1:length(goi_loci))){
  chrm = as.character(seqnames(goi_loci)[j])
  start = start(goi_loci)[j] - 20000
  end = end(goi_loci)[j] + 20000
  root = paste(sep = "_",
               chrm,
               paste0(as.integer(mean(c(start, end))/10^3), "K_maskFull"))
  all_best = list()
  pdf(paste0(root, ".pdf"))
  for(binned_name in names(all_binned)){
    key = binned_name
    keep = !(grepl("chrY", rownames(binned)) | grepl("chrM", rownames(binned)))
    binned = all_binned[[binned_name]]
    binned = binned[keep, keep]
    all_chrm = unique(sapply(strsplit(rownames(binned), " "), function(x)x[1]))
    names(all_chrm) = all_chrm
    chrm_sizes = sapply(all_chrm, function(x)grepl(paste0(x, " "), rownames(binned)) %>% sum)
    gene_name = goi_loci$gene_name[j]
    ucsc = paste0(as.character(seqnames(goi_loci))[j], ":", start(goi_loci)[j], "-", end(goi_loci)[j])
    all_best[[key]] = plot_region_interchromosomal(binned = binned, 
                                                   bin_size = bin_size,
                                                   main = paste0(key, ": ", gene_name, "\n", ucsc),
                                                   chrm = chrm, 
                                                   start = start, 
                                                   end = end, 
                                                   mask_n = 50000,#max(chrm_sizes), 
                                                   best_n = 1000)
  }
  save(all_best, file = paste0(root, ".save"))
  dev.off()
}
