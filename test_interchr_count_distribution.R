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

all_density_mats = list()


for(f in files){
  key = (((basename(f) %>% strsplit(split = "\\["))[[1]][2] %>% strsplit(split = "\\]"))[[1]][1] %>% strsplit(split = "_"))[[1]][1]
  print(f)
  load(f)
  all_chrm = unique(sapply(strsplit(rownames(binned), " "), function(x)x[1]))
  all_chrm = all_chrm[all_chrm != "chrY"]
  all_chrm = all_chrm[all_chrm != "chrM"]
  density_mat = matrix(0, ncol = length(all_chrm), nrow = length(all_chrm))
  colnames(density_mat) = all_chrm
  rownames(density_mat) = all_chrm
  chrm_mats = list()
  for(i in 1:length(all_chrm)){
    chrA = all_chrm[i]  
    
    
    # for(chrA in all_chrm){
    
    
    is_chrA  = grepl(paste0(chrA, " "), rownames(binned))
    # chr_binned = binned[is_chrA, is_chrA]
    for(j in i:length(all_chrm)){
      chrB = all_chrm[j]
      # for(chrB in all_chrm){
      print(paste(chrA, "->", chrB))
      is_chrB = grepl(paste0(chrB, " "), rownames(binned))
      chr_binned = binned[is_chrA, is_chrB]
      area = nrow(chr_binned) * ncol(chr_binned)
      total = sum(chr_binned)
      density_mat[chrA, chrB] = total / area
      density_mat[chrB, chrA] = total / area
    }
  }
  all_density_mats[[key]] = density_mat
}
# save(density_mat, file = "interchromosomal_interaction_density_heatmap.pdf")

save(all_density_mats, file = "all_interchromosomal_matrix.save")


colors = rgb(colorRamp(c("white", "gray", "blue", "purple", "red"))(0:50/50)/255)

pdf("interchromosomal_frequency.pdf")
hidden = sapply(names(all_density_mats), function(n){
  # nam = strsplit()
  x = all_density_mats[[n]]
  # heatmap.2(x, trace = "n")
  diag(x) = 0
  par(cex.main = .7)
  heatmap.2(x, trace = "n", main = n, col = colors)
})
dev.off()
