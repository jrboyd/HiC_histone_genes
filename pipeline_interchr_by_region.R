#summarizes interchromosomal interactions focused on a particular region
library(magrittr)
library(GenomicRanges)
source("functions_my_HiC_plots.R")

href_file = "hist_ref.save"
if(file.exists(href_file)){
  load(href_file)
}else{
  ref = parse_gtf("C:/Users/jrboyd/Downloads/gencode.v19.annotation.gtf", rownames_attrib = "gene_id", feature_type = "gene")
  is_hist = grepl("^HIST", ref$gene_name)
  hist_ref = ref[is_hist,]
  ref[ref$gene_name == "NPAT",]
  hist_ref = rbind(hist_ref, ref[ref$gene_name == "NPAT",], ref[ref$gene_name == "HINFP" ,])
  save(hist_ref, file = href_file)
}
loci_gap = 10^6
hist_pos = GRanges(seqnames = hist_ref$chrm, IRanges(hist_ref$start - loci_gap, hist_ref$end + loci_gap))
names(hist_pos) = hist_ref$gene_id
hist_loci = reduce(hist_pos)



files = dir(path = "bins_balanced_filtered2/", full.names = T)

plot_region_interchromosomal = function(binned, bin_size, main, chrm, start, end, mask_n, best_n){
  extract = extract_HiC_chrmRange(M = binned, chrm = chrm, start = start, end = end, bin_size = bin_size, by = "r")
  extract_means =  apply(extract, 2, mean)
  to_mask = genomic2matrix_index(start- bin_size * mask_n, end + bin_size * mask_n, bin_size)
  kept = setdiff(names(extract_means), paste(chrm, to_mask))
  extract_means = extract_means[kept]
  o = order(extract_means, decreasing = T)
  
  best =  head(extract_means[o], n = best_n)
  best = best[best > 0]
  if(length(best) == 0){
    plot(0:1, 0:1)
    text(.5,.5, "no interaction")
    return(NULL)
  }
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

all_binned = list()

for(j in (1:length(hist_loci))){
  chrm = as.character(seqnames(hist_loci)[j])
  start = start(hist_loci)[j] - 40001
  end = end(hist_loci)[j] + 40001
  print(paste(chrm, start, end))
  root = paste(sep = "_",
               chrm,
               paste0(as.integer(mean(c(start, end))/10^3), "K_maskIntra"))
  all_best = list()
  pdf(paste0(root, ".pdf"))
  for(f in files){
    key = basename(f) %>% sub(pattern = ".save", replacement = "")
    print(f)
    if(is.null(all_binned[[key]])){
      load(f)
      all_binned[[key]] = binned
    }else{
      binned = all_binned[[key]]
    }
    
    keep = !(grepl("chrY", rownames(binned)) | grepl("chrM", rownames(binned)))
    binned = binned[keep, keep]
    all_chrm = unique(sapply(strsplit(rownames(binned), " "), function(x)x[1]))
    names(all_chrm) = all_chrm
    chrm_sizes = sapply(all_chrm, function(x)grepl(paste0(x, " "), rownames(binned)) %>% sum)
    all_best[[key]] = plot_region_interchromosomal(binned = binned, 
                                                   bin_size = bin_size,
                                                   main = key,
                                                   chrm = chrm, 
                                                   start = start, 
                                                   end = end, 
                                                   mask_n = max(chrm_sizes), 
                                                   best_n = 1000)
  }
  save(all_best, file = paste0(root, ".save"))
  dev.off()
}
