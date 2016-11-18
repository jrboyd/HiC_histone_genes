library(Sushi)
library(GenomicRanges)
library(magrittr)
source("H:/R_workspace/jrb_R_scripts/parse_gtf.R")
source("make_bedpe.R")
if(file.exists("hist_ref.save")){
  load("hist_ref.save")
}else{
  ref = parse_gtf("C:/Users/jrboyd/Downloads/gencode.v24.annotation.gtf", rownames_attrib = "gene_id", feature_type = "gene")
  is_hist = grepl("^HIST", ref$gene_name)
  hist_ref = ref[is_hist,]
  ref[ref$gene_name == "NPAT",]
  hist_ref = rbind(hist_ref, ref[ref$gene_name == "NPAT",], ref[ref$gene_name == "HINFP" ,])
  save(hist_ref, file = "hist_ref.save")
}
files = c(
  "balanced_bins/Galaxy3-[hESC_balanced_bin_granges_40kb].data",
  "balanced_bins/Galaxy5-[MCF10A_balanced_bin_granges_40kb].data",
  "balanced_bins/Galaxy1-[MCF10A-AT1_balanced_bin_granges_40kb].data",
  "balanced_bins/Galaxy2-[MCF10A-CA1a_balanced_bin_granges_40kb].data",
  "balanced_bins/Galaxy4-[MCF7_balanced_bin_granges_40kb].data",
)
files = files[c(2,5)]
# files = files[c(1,3:4)]

all_binned = list()
for(f in files){
  key = (((basename(f) %>% strsplit(split = "\\["))[[1]][2] %>% strsplit(split = "\\]"))[[1]][1] %>% strsplit(split = "_"))[[1]][1]
  print(f)
  load(f)
  all_binned[[key]] = binned
  
}

loci_gap = 10^5
hist_pos = GRanges(seqnames = hist_ref$chrm, IRanges(hist_ref$start - loci_gap, hist_ref$end + loci_gap))
names(hist_pos) = hist_ref$gene_id
hist_loci = reduce(hist_pos)

#p_ext - how far plot is extended around loci
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
    
    k = grepl(roi.plotted$chrm, rownames(binned))
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
    phic = plotHic(log10(as.matrix(m*10000)+1), chrom = roi.plotted$chrm, chromstart = roi.plotted$start, chromend = roi.plotted$end, palette = colors, zrange = c(0,2.5), max_y = sqrt(nrow(m)^2/2))
    title(binned_name)
    axis(side = 1, at = roi.genes.pos, labels = roi.genes.sym, line = 1)
    labelgenome(roi$chrm,roi.plotted$start,roi.plotted$end,side=1,scipen=20,n=4,scale="Mb",edgeblankfraction=0.20,chromline=.5,scaleline=0.5, line = 3)
    addlegend(phic[[1]],palette=phic[[2]],title="score",side="right",bottominset=0.4,topinset=0,xoffset=-.035,labelside="left",width=0.025,title.offset=0.035)
  }
}

setwd("output_hmaps")
ext = 3*10^6
for(i in 1:length(all_binned)){
  for(j in 1:length(hist_loci)){
    root = paste(sep = "_",
          as.character(seqnames(hist_loci)[j]),
          paste0(as.integer(mean(c(start(hist_loci)[j], end(hist_loci)[j]))/10^3), "K"),
          names(all_binned)[i])
    print(root)
    # png(paste0(root, ".png"), width = 800, height = 600)
    # my_plotHic(names(all_binned)[i], loci_gr = hist_loci[j], hist_ref, p_ext = ext)
    # # chrm = as.character(seqnames(hist_loci)[j])
    # # start = start(hist_loci)[j]
    # # end = end(hist_loci)[j]
    # # my_plotHic_lines(binned_name = names(all_binned)[i], 
    # #                  chrm = chrm, start = start - ext, end = end + ext, 
    # #                  gene_ref = hist_ref, 
    # #                  min_dist = 4, 
    # #                  n_displayed = 200)
    # dev.off()
  }
  loci = hist_loci[11]
  end(loci) = end(hist_loci[12])
  png(paste0("NPAT-and-HINFP_", names(all_binned)[i], ".png"))
  my_plotHic(binned_name = names(all_binned)[i], loci_gr = loci, gene_ref = hist_ref, p_ext = ext)
  dev.off()
  # chrm = as.character(seqnames(loci))
  # start = start(loci)
  # end = end(loci)
  # my_plotHic_lines(binned_name = names(all_binned)[i],
  #                  chrm = chrm, start = start - ext, end = end + ext,
  #                  gene_ref = hist_ref,
  #                  min_dist = 10,
  #                  n_displayed = 200)
  # dev.off()
}
