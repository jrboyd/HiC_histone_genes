
# source("functions_moving_average.R")
my_plotHic = function(binned_name, chrm, start, end, gene_ref, applyLog = F, quantile_cap = 1, zrange = NULL){
  binned = all_binned[[binned_name]]
  
  # names(loci) = c("HIST2.1", "HIST2.2", "HIST2.3")
  # for(i in 1:length(loci_gr)){
    #range of interest
    roi = list(chrm = chrm, start = start, end = end)
    
    k = grepl(paste0(roi$chrm, " "), rownames(binned))
    m = binned[k,k]
    s = as.integer(roi$start / bin_size)
    e = as.integer(roi$end / bin_size) - 1
    m = m[s:e, s:e]
#rownames are starting coordinate so entry extends from rownames to rownames + bin_size
    rownames(m) = (s:e) * bin_size
    colnames(m) = rownames(m)
    if(max(m) ==0){
      plot(0:1, 0:1, type ="n")
      text(.5,.5, "no data to plot")
      return()
    }
    
    roi.gr = GRanges(seqnames = roi$chrm, IRanges(roi$start, roi$end))
    
    colors = function(n)rgb(colorRamp(c("white", "blue", "red"))((1:n-1)/n)/255)
    if(applyLog){
      mat = log10(as.matrix(m*10000)+1)
    }else{
      mat = as.matrix(m)
    }
    
    if(quantile_cap < 1){
      m_max = quantile(mat, quantile_cap)
      mat = ifelse(mat > m_max, m_max, mat)
    }
    
    
    if(is.null(zrange)) zrange = range(mat)
    phic = plotHic(mat, chrom = roi$chrm, chromstart = roi$start, chromend = roi$end, palette = colors, max_y = sqrt(nrow(m)^2/2), zrange = zrange)  
    
    title(binned_name)
    
    olaps = findOverlaps(query = roi.gr, subject = gene_ref)
    # roi.genes = names(gene_ref[subjectHits(olaps)])
    roi.gene_loci = gene_ref[subjectHits(olaps)]
    n_lines = 1
    for(i in 1:length(roi.gene_loci)){
      roi.gene_sym = roi.gene_loci$gene_name[i]
      axis(side = 1, at = c(start(roi.gene_loci)[i], end(roi.gene_loci)[i]), labels = rep("",2), line = i%%n_lines*2+0)
      axis(side = 1, at = mean(c(start(roi.gene_loci)[i], end(roi.gene_loci)[i])), labels = roi.gene_loci$gene_name[i], line = i%%n_lines*2+0, tick = F, las = 2, cex.axis = .2)
    }
    labelgenome(roi$chrm,roi$start,roi$end,side=1,scipen=20,n=4,scale="Mb",edgeblankfraction=0.20,chromline=.5,scaleline=0.5, line = 3, lwd = 3)
    axis(side = 1, at = 10^5 * 1:1000000, line = 3, cex = .5, labels = rep("", 1000000))
    addlegend(phic[[1]],palette=phic[[2]],title="score",side="right",bottominset=0.4,topinset=0,xoffset=-.035,labelside="left",width=0.025,title.offset=0.035)
  # }
}

my_plotBdg = function(binned_name, chrm, start, end, gene_ref){
  plot_start = start
  plot_end = end
  bedgraph_dir = paste0(wd, "/bedgraph_CTCF_FE/", binned_name)
  if(!dir.exists(bedgraph_dir))return()
  bedgraph_files = dir(bedgraph_dir, pattern = paste0(chrm, ".bdg"), full.names = T)
  if(length(bedgraph_files) == 0)return()
  bdg_i = 1
  is_multi = length(bedgraph_files) > 1
  colors = RColorBrewer::brewer.pal(max(length(bedgraph_files),3), "Set1")[1:length(bedgraph_files)]
  if(length(bedgraph_files) > 1) colors = paste0(colors, "99")
  for(bdg_f in bedgraph_files){
    bdg = read.table(bdg_f)
    head(bdg)
    colnames(bdg) = c("chrom", "start", "end", "value")
    bdg_gr = GRanges(bdg$chrom, IRanges(bdg$start, bdg$end), value = bdg$value)
    #sample every 100 bp
    end(bdg_gr) = end(bdg_gr) - 1
    s = as.integer(start / 100)+1
    e = as.integer(end / 100)
    xs = s:e*100
    bin_mids = xs
    bin_gr = GRanges(rep(chrm, length(xs)), IRanges(xs, xs))
    olaps = findOverlaps(bin_gr, bdg_gr)
    # bin_gr[queryHits(olaps)]
    xy_mat = cbind(xs, rep(0, length(xs)))
    rownames(xy_mat) = xs
    xy_mat[as.character(start(bin_gr)[queryHits(olaps)]), 2] = bdg_gr[subjectHits(olaps)]$value
    n = 100
    if(n > nrow(xy_mat)-1) n = nrow(xy_mat)-1
    size = as.integer((nrow(xy_mat)-1) / n)
    vmeans =  sapply(1:n, function(i){
      s = (i-1) * size + 1
      e = (i) * size
      e = ifelse(e > nrow(xy_mat), nrow(xy_mat), e)
      v = mean(xy_mat[s:e, 2])
    })
    poly_xs = c(rep(xs[length(xs)], 2), rep(xs[1], 2))
    poly_ys = c(vmeans[length(vmeans)], 0, 0, vmeans[1])
    for(i in 1:(length(vmeans)-1)){
      s = (i-1)*size + 1
      e = i * size
      poly_xs = c(poly_xs, xs[s])
      poly_xs = c(poly_xs, xs[e])
      poly_ys = c(poly_ys, vmeans[i])
      poly_ys = c(poly_ys, vmeans[i])
    }
    if(bdg_i == 1){
      plot(range(xs), c(0, max(vmeans)), type = "n", xaxs = "i", axes = F, xlab = "", ylab = "binned enrichment")
      legend("topright",  legend = basename(bedgraph_files) %>% sub(pattern = "_FE.+bdg", replacement = ""), fill = colors, bty = "n")
    }
    polygon(poly_xs, poly_ys, lty = 0, col = colors[bdg_i])
    bdg_i = bdg_i + 1
  }
  
  roi = list(chrm = chrm, start = start, end = end)
  roi.gr = GRanges(seqnames = roi$chrm, IRanges(roi$start, roi$end))
  olaps = findOverlaps(query = roi.gr, subject = gene_ref)
  # roi.genes = names(gene_ref[subjectHits(olaps)])
  roi.gene_loci = gene_ref[subjectHits(olaps)]
  n_lines = 1
  for(i in 1:length(roi.gene_loci)){
    roi.gene_sym = roi.gene_loci$gene_name[i]
    axis(side = 1, at = c(start(roi.gene_loci)[i], end(roi.gene_loci)[i]), labels = rep("",2), line = i%%n_lines*2+0)
    axis(side = 1, at = mean(c(start(roi.gene_loci)[i], end(roi.gene_loci)[i])), labels = roi.gene_loci$gene_name[i], line = i%%n_lines*2+0, tick = F, las = 2, cex.axis = .2)
  }
  labelgenome(roi$chrm,roi$start,roi$end,side=1,scipen=20,n=4,scale="Mb",edgeblankfraction=0.20,chromline=.5,scaleline=0.5, line = 3, lwd = 3)
  axis(side = 1, at = 10^5 * 1:1000000, line = 3, cex = .5, labels = rep("", 1000000))# roi$chrm,roi$start,roi$end,side=1,scipen=20,n=4,scale="100kb",edgeblankfraction=0.20,chromline=.5,scaleline=0.5, line = 3)
  axis(side = 2)
  return()
}

my_plotHic_lines = function(binned_name, chrm, start, end, gene_ref, bin_size = 40000, min_dist = 4, n_displayed = 100, quantile_cap = 1,zrange = NULL){
  bedpe = matrix2bedpe(M = all_binned[[binned_name]], chrm = chrm, start = start, end = end, min_dist = min_dist, bin_size = bin_size)
  n_displayed = min(n_displayed, nrow(bedpe))
  toplot = bedpe
  # MIN = quantile(toplot$score, min_quant)
  o = order(toplot$score, decreasing = T)
  if(quantile_cap < 1){
    score_max = quantile(toplot$score, quantile_cap)
    toplot$score = ifelse(toplot$score > score_max, score_max, toplot$score)
    
  }
  
  MIN = toplot$score[o][n_displayed]
  toplot = subset(toplot, score >= MIN)
  # toplot = bedpe
  toplot$distance = abs(toplot$start2 - toplot$start1)
  toplot$weight = toplot$distance / (max(min_dist,1) * bin_size * 2)+.1
  toplot$weight = ifelse(toplot$weight > 1, 1, toplot$weight)
  # toplot = subset(toplot, distance > bin_size * 4)
  if(!is.null(zrange)) toplot$score = ifelse(toplot$score > max(zrange), max(zrange), toplot$score)
  plotBedpe(bedpedata = toplot[,1:8], chrom = chrm, chromstart = start, chromend = end, heights = toplot$score, plottype = "loops", xaxt = "n")
  # labelgenome(chrm, start, end, n=3,scale="Mb")
  # legend("topright",inset =0.01,legend=c("K562","HeLa","GM12878"),
  #        col=SushiColors(3)(3),pch=19,bty="n",text.font=2)
  axis(side=2,las=2,tcl=.2)
  title(binned_name)
  roi = list(chrm = chrm, start = start, end = end)
  roi.gr = GRanges(seqnames = chrm, IRanges(start, end))
  olaps = findOverlaps(query = roi.gr, subject = gene_ref)
  # roi.genes = names(gene_ref[subjectHits(olaps)])
  roi.gene_loci = gene_ref[subjectHits(olaps)]
  n_lines = 1
  for(i in 1:length(roi.gene_loci)){
    roi.gene_sym = roi.gene_loci$gene_name[i]
    axis(side = 1, at = c(start(roi.gene_loci)[i], end(roi.gene_loci)[i]), labels = rep("",2), line = i%%n_lines*2+0)
    axis(side = 1, at = mean(c(start(roi.gene_loci)[i], end(roi.gene_loci)[i])), labels = roi.gene_loci$gene_name[i], line = i%%n_lines*2+0, tick = F, las = 2, cex.axis = .2)
  }
  labelgenome(roi$chrm,roi$start,roi$end,side=1,scipen=20,n=4,scale="Mb",edgeblankfraction=0.20,chromline=.5,scaleline=0.5, line = 3, lwd = 3)
  axis(side = 1, at = 10^5 * 1:1000000, line = 3, cex = .5, labels = rep("", 1000000))# roi$chrm,roi$start,roi$end,side=1,scipen=20,n=4,scale="100kb",edgeblankfraction=0.20,chromline=.5,scaleline=0.5, line = 3)
  
  roi.gene_pos = rowMeans(cbind(start(roi.gene_loci), end(roi.gene_loci)))
  # axis(side = 1, at = roi.gene_pos, labels = roi.gene_sym, line = 1)
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
  s = as.integer(roi$start / bin_size) 
  e = as.integer(roi$end / bin_size) -1
  m = m[s:e, s:e]
  rownames(m) = (s:e) * bin_size
  colnames(m) = rownames(m)
  
  #sort be decreasing strength
  o = order(m@x, decreasing = T)
  mat = data.frame(x = m@x, i = m@i+s, j = m@j+s)
  mat = mat[o,]
  
  #limit to over minimum distance
  mat$dist = mat$j - mat$i
  keep = mat$dist >= min_dist
  mat = mat[keep,]
  
  #generate bedpe
  len = nrow(mat)
  bedpe = data.frame(chrom1 = rep(chrm, len), start1 = (mat$i-0) * bin_size, end1 = (mat$i+1) * bin_size, 
                     chrom2 = rep(chrm, len), start2 = (mat$j-0) * bin_size, end2 = (mat$j+1) * bin_size,
                     name = rep(NA, len), score = mat$x,
                     strand1 = rep(".", len), strand2 = rep(".", len),
                     samplenumber = rep(samplenumber, len))
  
  # #old insane method for generating bedpe
  # # bedpe = character(nrow(m) * (ncol(m)+1)/2)
  # size = nrow(m) - min_dist
  # len = size * (size+1)/2
  # bedpe = data.frame(chrom1 = rep(chrm, len), start1 = integer(len), end1 = integer(len), 
  #                    chrom2 = rep(chrm, len), start2 = integer(len), end2 = integer(len), 
  #                    name = rep(NA, len), score = numeric(len), 
  #                    strand1 = rep(".", len), strand2 = rep(".", len), 
  #                    samplenumber = rep(samplenumber, len)) 
  # a = 1
  # for(i in 1:(nrow(m)-min_dist)){
  #   for(j in (i + min_dist):ncol(m)){
  #     # bedpe[a] = paste(sep = "\t",
  #     #                  chrm, as.integer(rownames(m)[i]) - bin_size, rownames(m)[i],
  #     #                  chrm, as.integer(colnames(m)[j]) - bin_size, colnames(m)[j],
  #     #                  NA, m[i,j], ".", ".", 1)
  #     bedpe[a,]$start1 = as.integer(rownames(m)[i]) - bin_size
  #     bedpe[a,]$end1 = as.integer(rownames(m)[i])
  #     bedpe[a,]$start2 = as.integer(colnames(m)[j]) - bin_size
  #     bedpe[a,]$end2 = as.integer(colnames(m)[j])
  #     bedpe[a,]$score = m[i,j]
  #     a = a + 1
  #   }
  #   
  # }
  return(bedpe)
}

