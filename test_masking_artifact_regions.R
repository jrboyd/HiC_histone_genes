#masking regions <.1 quantile for regions interacted with seems appropriate
library(magrittr)
#characterize low mapping regions that are being messed up in bin balancing
# load("balanced_bins/Galaxy3-[AT1_TotUnbalMatrix_hg19_40kb].data")
files = dir(path = "bins_balanced/", full.names = T)
# files = files[c(2,5)]
# # files = files[c(1,3:4)]

pdf("count_hists.pdf")
all_binned = list()
for(f in files){
  key = (((basename(f) %>% strsplit(split = "\\["))[[1]][2] %>% strsplit(split = "\\]"))[[1]][1] %>% strsplit(split = "_"))[[1]][1]
  print(f)
  load(f)
  is_M = grepl("chrM", rownames(binned))
  binned = binned[!is_M, !is_M]
  cSums = colSums(binned)
  cCounts = colSums(binned > 0)
  
  MIN = quantile(cCounts, .15)
  k = cCounts >= MIN
  
  hist((cCounts[k]), xlim = c(0,3000), breaks = max(cCounts)+1, lty = 0, col = "black")
  lines(hist((cCounts[!k]), breaks = MIN+1, plot = F), col = "red", lty = 0)  
  
}
dev.off()
