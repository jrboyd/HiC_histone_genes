library(data.table)
library(GenomicRanges)
library(Matrix)
#proper table write format for tadtools
my_write.table = function(dat, file){
  write.table(dat, file = file, sep ="\t", row.names = F, col.names = F, quote = F)
}
#load and save example data, known to work
# ex_tadtool_matrix_file = "~/tadtool/examples/chr12_20-35Mb.matrix.txt"
# ex_tadtool_bed_file = "~/tadtool/examples/chr12_20-35Mb_regions.bed"
# ex_tadtool_matrix_mat = read.table(ex_tadtool_matrix_file)
# ex_tadtool_bed_mat = read.table(ex_tadtool_bed_file)
# my_write.table(ex_tadtool_matrix_mat,"tads_input/ex_matrix.txt")
# my_write.table(ex_tadtool_bed_mat, "tads_input/ex_regions.bed")

convert_HiCPro_OUT_to_WUbrowser_IN = function(hicp_matrix_file, hicp_bed_file, cell, binSize, chr, out_dir = "~/R/HiC_histone_genes/wubrowser2_format/"){
  prefix = paste0(cell, "_", binSize, "bin")
  
  #load real data output by HiC-Pro and convert to appropriate formats
  hicp_matrix_dt = fread(hicp_matrix_file)
  colnames(hicp_matrix_dt) = c("r1", "r2", "val")
  hicp_bed_mat = read.table(hicp_bed_file)
  colnames(hicp_bed_mat) = c("seqnames", "start", "end", "id")
  # hicp_bed_gr = GRanges(hicp_bed_mat)
  chr_bed_mat = subset(hicp_bed_mat, seqnames == chr)
  chr_matrix_dt = hicp_matrix_dt[r1  %in% chr_bed_mat$id & r2 %in% chr_bed_mat$id]
  
  a = chr_bed_mat[as.character(chr_matrix_dt$r1),]
  colnames(a) = paste0("r1_", colnames(a))
  a = as.data.table(a)
  b = chr_bed_mat[as.character(chr_matrix_dt$r2),]
  colnames(b) = paste0("r2_", colnames(b))
  b = as.data.table(b)
  chr_matrix_dt = cbind(chr_matrix_dt, a, b)
  # tmp = head(chr_matrix_dt)
  joined = chr_matrix_dt[, .(paste0(r1_seqnames, ":", r1_start, "-", r1_end), paste0(r2_seqnames, ":", r2_start, "-", r2_end), val)]  
  
  dir.create(out_dir, showWarnings = F)
  start_dir = getwd()
  setwd(out_dir)
  # df = as.data.frame(topL)
  my_write.table(joined, file = paste0(prefix, "_matrix.txt"))
  setwd(start_dir)
  return(joined)
}
#pooled samples
setwd("~/HiC-Pro/outputs/MCF10A_pooled/hic_results/matrix/MCF10A_20split")
M = convert_HiCPro_OUT_to_WUbrowser_IN(hicp_matrix_file = "iced/40000/MCF10A_20split_40000_iced.matrix", 
                                       hicp_bed_file = "raw/40000/MCF10A_20split_40000_abs.bed",
                                       cell = "MCF10A_pooled",
                                       binSize = 40000, 
                                       chr = "chr6")
# , 
# chr = "chr6", start = 23*10^6, end = 31*10^6)
for(cl in c("MCF10AT1", "MCF10CA1a")){
  for(n in c("pooled")){
    print(paste(cl, n))
    setwd(paste0("~/HiC-Pro/outputs/", cl, "_pooled/hic_results/matrix/", n))
    convert_HiCPro_OUT_to_WUbrowser_IN(hicp_matrix_file = paste0("iced/40000/", n, "_40000_iced.matrix"),
                                       hicp_bed_file = paste0("raw/40000/", n, "_40000_abs.bed"),
                                       cell = paste0(cl, "_", n),
                                       binSize = 40000,
                                       chr = "chr6")
  }
}
