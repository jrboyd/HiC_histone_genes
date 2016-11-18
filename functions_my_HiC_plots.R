
my_plotHic = function(binned_name, loci_gr, gene_ref, p_ext = 5*10^5){
  binned = all_binned[[binned_name]]
  
  # names(loci) = c("HIST2.1", "HIST2.2", "HIST2.3")
  for(i in 1:length(loci_gr)){
    #range of interest
    roi = list(chrm = as.character(seqnames(loci_gr))[i], start = start(loci_gr)[i], end = end(loci_gr)[i])
    roi.plotted = roi #range to be plotted
    
    roi.plotted$start = max(roi.plotted$start - p_ext, 0)
    roi.plotted$end = roi.plotted$end + p_ext
    # roi.plotted$start = 25000000
    # roi.plotted$end = 29000000
    
    k = grepl(paste0(roi.plotted$chrm, " "), rownames(binned))
    m = binned[k,k]
    s = as.integer(roi.plotted$start / bin_size)
    e = as.integer(roi.plotted$end / bin_size) + 1
    m = m[s:e, s:e]
    rownames(m) = (s:e - 1) * bin_size
    colnames(m) = rownames(m)
    
    roi.plotted.gr = GRanges(seqnames = roi.plotted$chrm, IRanges(roi.plotted$start, roi.plotted$end))
    olaps = findOverlaps(query = roi.plotted.gr, subject = hist_pos)
    roi.genes = names(hist_pos[subjectHits(olaps)])
    roi.genes.sym = gene_ref[roi.genes,]$gene_name
    roi.genes.pos = rowMeans(gene_ref[roi.genes,4:5])
    
    colors = function(n)rgb(colorRamp(c("white", "blue", "red"))((1:n-1)/n)/255)
    phic = plotHic(log10(as.matrix(m*10000)+1), chrom = roi.plotted$chrm, chromstart = roi.plotted$start, chromend = roi.plotted$end, palette = colors, max_y = sqrt(nrow(m)^2/2))
    title(binned_name)
    axis(side = 1, at = roi.genes.pos, labels = roi.genes.sym, line = 1)
    labelgenome(roi$chrm,roi.plotted$start,roi.plotted$end,side=1,scipen=20,n=4,scale="Mb",edgeblankfraction=0.20,chromline=.5,scaleline=0.5, line = 3)
    addlegend(phic[[1]],palette=phic[[2]],title="score",side="right",bottominset=0.4,topinset=0,xoffset=-.035,labelside="left",width=0.025,title.offset=0.035)
  }
}

my_plotHic_lines = function(binned_name, chrm, start, end, gene_ref, min_dist = 4, n_displayed = 1000){
  bedpe = matrix2bedpe(all_binned[[binned_name]], chrm = chrm, start = start, end = end, min_dist = min_dist)
  toplot = bedpe
  # MIN = quantile(toplot$score, min_quant)
  o = order(toplot$score, decreasing = T)
  MIN = toplot$score[o][n_displayed]
  toplot = subset(toplot, score > MIN)
  # toplot = bedpe
  toplot$distance = abs(toplot$start2 - toplot$start1)
  toplot = subset(toplot, distance > bin_size * 4)
  plotBedpe(bedpedata = toplot, chrom = chrm, chromstart = start, chromend = end, heights = toplot$score, plottype = "loops", lwdby = toplot$distance)
  labelgenome(chrm, start, end, n=3,scale="Mb")
  # legend("topright",inset =0.01,legend=c("K562","HeLa","GM12878"),
  #        col=SushiColors(3)(3),pch=19,bty="n",text.font=2)
  axis(side=2,las=2,tcl=.2)
  title(binned_name)
  roi.gr = GRanges(seqnames = chrm, IRanges(start, end))
  olaps = findOverlaps(query = roi.gr, subject = hist_pos)
  roi.genes = names(hist_pos[subjectHits(olaps)])
  roi.genes.sym = gene_ref[roi.genes,]$gene_name
  roi.genes.pos = rowMeans(gene_ref[roi.genes,4:5])
  axis(side = 1, at = roi.genes.pos, labels = roi.genes.sym, line = 1)
}

genomic2matrix_index = function(start, end, bin_size){
  s = as.integer(start / bin_size) + 1
  e = as.integer((end-1) / bin_size) + 1
  return(s:e)
}

extract_HiC_chrmRange = function(M, chrm, start, end, bin_size, by = "rc"){
  extract_HiC_chrm(M, chrm, by) %>% extract_HiC_range(start, end, bin_size, by)
}

#M is assumed to be a square matrix symetrically representing genomic locations with rownames like "chr1 1, chr1 2, etc."
extract_HiC_chrm = function(M, chrm, by = "rc"){
  do_row = grepl("r", by) 
  do_col = grepl("c", by)
  if(! (do_row || do_col)) stop("by must contain r and/or c to indicate row and/or column")
  k = grepl(paste0(chrm, " "), rownames(M))
  kr = T
  kc = T
  if(do_row) kr = k
  if(do_col) kc = k
  return(M[kr,kc])
}

#M is assumed to be a square matrix symetrically representing genomic locations with rownames like "chr1 1, chr1 2, etc."
#M must represent a single chr along the dimension indicated by by (run extract_HiC_chrm first)
extract_HiC_range = function(M, start, end, bin_size, by = "rc"){
  do_row = grepl("r", by) 
  do_col = grepl("c", by)
  if(! (do_row || do_col)) stop("by must contain r and/or c to indicate row and/or column")
  if(do_row){
    if(sapply(strsplit(rownames(M), " "), function(x)x[1]) %>% unique %>% length > 1)
      stop("rows of M contains multiple chrm, run extract_HiC_chrm")
  }
  if(do_col){
    if(sapply(strsplit(colnames(M), " "), function(x)x[1]) %>% unique %>% length > 1)
      stop("columns of M contains multiple chrm, run extract_HiC_chrm")
  }
  s = as.integer(start / bin_size) + 1
  e = as.integer((end-1) / bin_size) + 1
  kr = T
  kc = T
  if(do_row) kr = s:e
  if(do_col) kc = s:e
  m = M[kr, kc]
  options(scipen=999)
  if(do_row) rownames(m) = (kr - 1) * bin_size
  if(do_col) colnames(m) = (kc - 1) * bin_size
  options(scipen=0)
  return(m)
}

matrix2bedpe = function(M, chrm, start, end, bin_size = 40000, samplenumber = 1, min_dist = 1){
  options(scipen=999)
  roi = list(chrm = chrm, start = start, end = end)
  
  k = grepl(paste0(roi$chrm, " "), rownames(M))
  m = M[k,k]
  s = as.integer(roi$start / bin_size) + 2
  e = as.integer(roi$end / bin_size) + 1
  m = m[s:e, s:e]
  rownames(m) = (s:e - 1) * bin_size
  colnames(m) = rownames(m)
  
  # bedpe = character(nrow(m) * (ncol(m)+1)/2)
  size = nrow(m) - min_dist
  len = size * (size+1)/2
  bedpe = data.frame(chrom1 = rep(chrm, len), start1 = integer(len), end1 = integer(len), 
                     chrom2 = rep(chrm, len), start2 = integer(len), end2 = integer(len), 
                     name = rep(NA, len), score = numeric(len), 
                     strand1 = rep(".", len), strand2 = rep(".", len), 
                     samplenumber = rep(samplenumber, len)) 
  a = 1
  for(i in 1:(nrow(m)-min_dist)){
    for(j in (i + min_dist):ncol(m)){
      # bedpe[a] = paste(sep = "\t",
      #                  chrm, as.integer(rownames(m)[i]) - bin_size, rownames(m)[i],
      #                  chrm, as.integer(colnames(m)[j]) - bin_size, colnames(m)[j],
      #                  NA, m[i,j], ".", ".", 1)
      bedpe[a,]$start1 = as.integer(rownames(m)[i]) - bin_size
      bedpe[a,]$end1 = as.integer(rownames(m)[i])
      bedpe[a,]$start2 = as.integer(colnames(m)[j]) - bin_size
      bedpe[a,]$end2 = as.integer(colnames(m)[j])
      bedpe[a,]$score = m[i,j]
      a = a + 1
    }
    
  }
  return(bedpe)
}

