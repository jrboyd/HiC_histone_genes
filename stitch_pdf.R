library(png)
wd = getwd()
setwd("output_AF_hg38_integrated_filtered2_ctcf_summary/")
png_files = dir(pattern = ".png")
roots = sapply(strsplit(png_files, "_"), function(x)paste(x[1:2], collapse = "_"))
roots_uniq = unique(roots)
for(root in roots_uniq){
  k = root == roots
  to_stitch = png_files[k]
  to_stitch = to_stitch[c(3,1,2)]
  pdf(paste0(root, ".pdf"), width = 8, height = 4)
  par(mai = rep(.1,4))
  layout(matrix(1:3, nrow = 1, byrow = T))
  for(png_f in to_stitch){
    cl = strsplit(png_f, "[_\\.]")[[1]][3]
    plot(0, xlim = 0:1, ylim = 0:1, type = "n", xlab = "", ylab = "", axes = F)
    title(cl)
    rasterImage(readPNG(png_f, native = FALSE), xleft = 0, xright = 1, ytop = 1, ybottom = 0,
                interpolate = FALSE) 
  }
  dev.off()
}