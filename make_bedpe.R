# data(package = "Sushi", "Sushi_5C.bedpe")
# head(Sushi_5C.bedpe)

# M = all_binned[[1]]
# chrm = "chr6"
# start = 25000000
# end = 26000000

matrix2bedpe = function(M, chrm, start, end, bin_size = 40000, samplenumber = 1, min_dist = 1){
  options(scipen=999)
  roi = list(chrm = chrm, start = start, end = end)
  
  k = grepl(roi$chrm, rownames(M))
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

# chrm = "chr6"
# start = 25*10^6
# end = 29*10^6
# gene_ref = hist_ref


# mtext("Z-score",side=2,line=1.75,cex=.75,font=2)
