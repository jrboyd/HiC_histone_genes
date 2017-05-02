library(Sushi)
library(GenomicRanges)
library(magrittr)
source("parse_gtf.R")
source("functions_my_HiC_plots.R")
if(file.exists("hist_ref.save")){
  load("hist_ref.save")
}else{
  ref = parse_gtf("gencode.v24.annotation.gtf", rownames_attrib = "gene_id", feature_type = "gene")
  is_hist = grepl("^HIST", ref$gene_name)
  hist_ref = ref[is_hist,]
  ref[ref$gene_name == "NPAT",]
  hist_ref = rbind(hist_ref, ref[ref$gene_name == "NPAT",], ref[ref$gene_name == "HINFP" ,])
  save(hist_ref, file = "hist_ref.save")
}
loci_gap = 10*10^6
hist_pos = GRanges(seqnames = hist_ref$chrm, IRanges(hist_ref$start - loci_gap, hist_ref$end + loci_gap))
names(hist_pos) = hist_ref$gene_id
hist_loci = reduce(hist_pos)


files = dir("bins_balanced_filtered2/", full.names = T)
# files = files[-1]
#   c(
#   "bins/MCF10A-AT1.save" #,
#   # "balanced_bins/Galaxy5-[MCF10A_balanced_bin_granges_40kb].data",
#   # "balanced_bins/Galaxy1-[MCF10A-AT1_balanced_bin_granges_40kb].data",
#   # "balanced_bins/Galaxy2-[MCF10A-CA1a_balanced_bin_granges_40kb].data",
#   # "balanced_bins/Galaxy4-[MCF7_balanced_bin_granges_40kb].data"
# )
# files = files[c(3)]
# files = files[c(1,3:4)]

all_binned = list()
for(f in files){
  key = basename(f) %>% sub(pattern = ".save", replacement = "")
  print(f)
  load(f)
  all_binned[[key]] = binned
}

my_plots = function(binned_name, loci_gr, gene_ref, ext){
  
  layout(1:3, heights = c(3,2,1))
  chrm = as.character(seqnames(loci_gr))
  start = start(loci_gr)
  end = end(loci_gr)
  
  my_plotHic(binned_name = binned_name, 
             loci_gr = loci_gr, 
             gene_ref = gene_ref, 
             p_ext = ext)
  my_plotHic_lines(binned_name = binned_name,
                   chrm = chrm, start = start - ext, end = end + ext,
                   gene_ref = gene_ref,
                   min_dist = 4,
                   n_displayed = 200)
  #do bedgraph plot
  plot_start = as.integer((start - ext) / bin_size) * bin_size
  plot_end = as.integer((end + ext + bin_size) / bin_size) * bin_size
  bedgraph_dir = paste0(wd, "/bedgraph_CTCF_FE/", names(all_binned)[i])
  if(!dir.exists(bedgraph_dir))return()
  bedgraph_files = dir(bedgraph_dir, pattern = paste0(chrm, ".bdg"), full.names = T)
  if(length(bedgraph_files) == 0)return()
  bdg_i = 1
  is_multi = length(bedgraph_files) > 1
  colors = RColorBrewer::brewer.pal(max(length(bedgraph_files),3), "Set1")[1:length(bedgraph_files)]
  for(bdg_f in bedgraph_files){
    bdg = read.table(bdg_f)
    head(bdg)
    colnames(bdg) = c("chrom", "start", "end", "value")
    bdg = bdg[bdg$value > 5,]
    # bdg$value = ifelse(bdg$value
    plotBedgraph(bdg, chrom = chrm, chromstart = plot_start, chromend = plot_end, 
                 overlay = bdg_i > 1, transparency = ifelse(is_multi, .6, 0), rescaleoverlay = is_multi, color = colors[bdg_i])
    if(bdg_i == 1){
      legend("topright",  legend = basename(bedgraph_files) %>% sub(pattern = "_FE.+bdg", replacement = ""), fill = colors, bty = "n")
      axis(1)
      axis(2)
    }
    bdg_i = bdg_i + 1
  }
  
}

wd = getwd()
outdir = normalizePath("output_hg38_integrated_filtered2")
dir.create(outdir)
setwd(outdir)

for(i in 1:length(all_binned)){
  for(j in 1:length(hist_loci)){
    root = paste(sep = "_",
                 as.character(seqnames(hist_loci)[j]),
                 paste0(as.integer(mean(c(start(hist_loci)[j], end(hist_loci)[j]))/10^3), "K"),
                 names(all_binned)[i])
    print(root)
    png(paste0(root, ".png"), width = 600, height = 1000)
    my_plots(binned_name = names(all_binned)[i], 
             loci_gr = hist_loci[j], 
             gene_ref = hist_ref, 
             ext = 2*10^5)
    dev.off()
  }
  print("NPAT and HINFP")
  loci = hist_loci[length(hist_loci)-1]
  end(loci) = end(hist_loci[length(hist_loci)])
  chrm = as.character(seqnames(loci))
  start = start(loci)
  end = end(loci)
  
  png(paste0("NPAT-and-HINFP_", names(all_binned)[i], ".png"), width = 800, height = 1800)
  my_plots(binned_name = names(all_binned)[i], 
           loci_gr = loci, 
           gene_ref = hist_ref, 
           ext = 2*10^5)
  dev.off()
}
setwd('..')


