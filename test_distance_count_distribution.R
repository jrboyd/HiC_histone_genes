#masking regions <.1 quantile for regions interacted with seems appropriate
library(magrittr)
#characterize low mapping regions that are being messed up in bin balancing
# load("balanced_bins/Galaxy3-[AT1_TotUnbalMatrix_hg19_40kb].data")
files = c(
  "balanced_bins/Galaxy4-[hESC_TotUnbalMatrix_hg19_40kb].data",
  "balanced_bins/Galaxy2-[10A_TotUnbalMatrix_hg19_40kb].data",
  "balanced_bins/Galaxy3-[AT1_TotUnbalMatrix_hg19_40kb].data",
  "balanced_bins/Galaxy1-[CA1a_TotUnbalMatrix_hg19_40kb].data",
  "balanced_bins/Galaxy5-[MCF7_TotUnbalMatrix_hg19_40kb].data"
)
# files = files[c(2,5)]
# # files = files[c(1,3:4)]
randi = function(n_rand, from){
  order(runif(from))[1:n_rand]
}

all_mats = list()


pdf("mat_profs.pdf")
for(f in files[1]){
  key = (((basename(f) %>% strsplit(split = "\\["))[[1]][2] %>% strsplit(split = "\\]"))[[1]][1] %>% strsplit(split = "_"))[[1]][1]
  print(f)
  load(f)
  all_chrm = unique(sapply(strsplit(rownames(binned), " "), function(x)x[1]))
  chrm_mats = list()
  for(chr in all_chrm){
    
    if(chr == "chrY" || chr == "chrM") next
    print(chr)
    is_chr  = grepl(paste0(chr, " "), rownames(binned))
    chr_binned = binned[is_chr, is_chr]
    MAX_DIAG = 100
    MAX_DIAG = min(MAX_DIAG, nrow(chr_binned))
    doit = function(){
      pb <- txtProgressBar(style = 3)
      setTxtProgressBar(pb, 0)
      mat = sapply(0:MAX_DIAG, function(i){#for each distance, starting at 0 (diagonal) then displacing by one, sample 300 interaction counts
        sample = randi(300, nrow(chr_binned)-i)
        # val = chr_binned[sample, sample+i]
        # val = sapply(randi(300, nrow(chr_binned)-i), function(x)chr_binned[x,x+i])
        if(i > 0){
          val = quantile(diag(chr_binned[,-1:-i]), c(.05,.25,.5,.75,.95))
        }else{
          val = quantile(diag(chr_binned), c(.05,.25,.5,.75,.95))
        }
        
        setTxtProgressBar(pb, i/MAX_DIAG)
        return(val)
      })
      close(pb)
      return(mat)
    }
    mat = doit()
    chrm_mats[[chr]] = mat
    plot(mat[5,], type = "l", xlim = c(0,50), log = "y")
    lines(mat[4,])
    lines(mat[3,])
    lines(mat[2,])
    title(paste(f, chr, sep = "\n"))
  }
  all_mats[[f]] = chrm_mats
}
remove(binned)
dev.off()

save(all_mats, file = "all_intrachromosomal_quantiles.save")

qi = 4
chrm = "chr6"
plot_intra_quantiles = function(chrm, qi = 4){
  q_name = rownames(all_mats$`balanced_bins/Galaxy2-[10A_TotUnbalMatrix_hg19_40kb].data`$chr1)[qi]
  vals = sapply(all_mats, function(x)x[[chrm]][qi,])
  colnames(vals) = sapply(strsplit(colnames(vals), "[\\[_]"), function(x)x[3])
  norm_vals = apply(vals, 2, function(x)x/x[1])
  plot(0, xlim = c(0,10), ylim = c(0,1), type= "n", ylab = "relative frequency", xlab = "bin distance")
  colors = RColorBrewer::brewer.pal(ncol(vals), "Set1")
  for(i in 1:ncol(norm_vals)){
    lines(x = 1:nrow(norm_vals)-1, norm_vals[,i], col = colors[i], lwd = 3)
  }
  legend(x = "topright",legend = colnames(vals), fill = colors)
  title(paste(q_name, "quantile of interactions intra", chrm))
}

plot_intra_quantiles("chr1")
