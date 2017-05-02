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

convert_HiCPro_OUT_to_Tadtool_IN = function(hicp_matrix_file, hicp_bed_file, cell, binSize, chr, start, end, out_dir = "~/R/HiC_histone_genes/tads_input"){
  prefix = paste0(cell, "_", binSize, "bin")
  
  #load real data output by HiC-Pro and convert to appropriate formats
  hicp_matrix_dt = fread(hicp_matrix_file)
  hicp_bed_mat = read.table(hicp_bed_file)
  colnames(hicp_bed_mat) = c("seqnames", "start", "end", "id")
  hicp_bed_gr = GRanges(hicp_bed_mat)
  
  #define target region of interest
  target_regions = GRanges(chr, IRanges(start+1, end-1))
  #intersect hicp bed with target region
  tadtool_bed_df = as.data.frame(hicp_bed_gr[queryHits(findOverlaps(hicp_bed_gr, target_regions))])
  #matrix ids to subset
  target_ids = tadtool_bed_df$id
  #tweak format to match tadtool example
  tadtool_bed_df[,4:6] = NULL
  tadtool_bed_df[,2] = tadtool_bed_df[,2] + 1
  options(scipen = 999)
  
  dir.create(out_dir, showWarnings = F)
  start_dir = getwd()
  setwd(out_dir)
  my_write.table(tadtool_bed_df, file = paste0(paste(prefix, chr, start, end, sep = "_"), "_regions.bed"))
  #select entries in matrix with V1 or V2 as target_id
  #data.table is awesome for this
  target_matrix_dt = hicp_matrix_dt[V1  %in% target_ids & V2 %in% target_ids]
  target_matrix_dt[, c("V1", "V2") := .(V1 - min(target_ids) , V2 - min(target_ids) )]
  diag = target_matrix_dt[V1 == V2]
  main = target_matrix_dt[V1 != V2]
  target_matrix_dt = rbind(diag, main, main[, .(V1 = V2, V2 = V1, V3 = V3)])
  # target_matrix_dt[, .(V1 = V1 - min(target_ids) , V2 = V2 - min(target_ids) )]
  M = new("dgTMatrix",
          i = target_matrix_dt$V1,
          j = target_matrix_dt$V2, 
          x=target_matrix_dt$V3, 
          Dim= c(length(target_ids), length(target_ids)))
  
  rowSums(M)
  my_write.table(as.matrix(M), file = paste0(paste(prefix, chr, start, end, sep = "_"), "_matrix.txt"))
  setwd(start_dir)
  return(M)
}
start = 23*10^6
end = 31*10^6
start = 0
end = 170805979
#pooled samples
setwd("~/HiC-Pro/outputs/MCF10A_pooled/hic_results/matrix/MCF10A_20split")
M = convert_HiCPro_OUT_to_Tadtool_IN(hicp_matrix_file = "iced/40000/MCF10A_20split_40000_iced.matrix", 
                                 hicp_bed_file = "raw/40000/MCF10A_20split_40000_abs.bed",
                                 cell = "MCF10A_pooled",
                                 binSize = 40000, 
                                 chr = "chr6", start = start, end = end)
for(cl in c("MCF10AT1", "MCF10CA1a")){
  for(n in c("pooled")){
    print(paste(cl, n))
    setwd(paste0("~/HiC-Pro/outputs/", cl, "_pooled/hic_results/matrix/", n))
    convert_HiCPro_OUT_to_Tadtool_IN(hicp_matrix_file = paste0("iced/40000/", n, "_40000_iced.matrix"),
                                     hicp_bed_file = paste0("raw/40000/", n, "_40000_abs.bed"),
                                     cell = paste0(cl, "_", n),
                                     binSize = 40000,
                                     chr = "chr6", start = start, end = end)
  }
}
#rep samples
for(cl in c("MCF10A", "MCF10AT1", "MCF10CA1a")){
  for(n in c("rep1", "rep2")){
    print(paste(cl, n))
    setwd(paste0("~/HiC-Pro/outputs/", cl, "/hic_results/matrix/", n))
    convert_HiCPro_OUT_to_Tadtool_IN(hicp_matrix_file = paste0("iced/40000/", n, "_40000_iced.matrix"), 
                                     hicp_bed_file = paste0("raw/40000/", n, "_40000_abs.bed"),
                                     cell = paste0(cl, "_", n),
                                     binSize = 40000, 
                                     chr = "chr6", start = 23*10^6, end = 31*10^6)
  }
}
