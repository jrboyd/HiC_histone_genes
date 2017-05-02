library(Sushi)
library(GenomicRanges)
library(magrittr)
source("parse_gtf.R")
source("functions_my_HiC_plots.R")

load("hg38ref.save")
if(F){
  ref = parse_gtf("gencode.v24.annotation.gtf", rownames_attrib = "gene_id", feature_type = "gene", additional_attrib = "gene_type")
  ref = ref[ref$gene_type == "protein_coding",]
  uniq_genes = unique(ref$gene_name)
  
  uniq_ref = ref[0,]
  for(g in uniq_genes){
    toadd = subset(ref, gene_name == g)
    if(length(unique(toadd$chrm)) > 1){
      print(g)
      next
    }
    uniq_ref = rbind(uniq_ref, subset(ref, gene_name == g))
  }
  ref = uniq_ref
  save(ref, file = "hg38ref.save")
}
# goi = c("ZEB1", "PLAUR", "IL6", "CDH1", "FN1")
# goi = "FN1"
# goi_ref = ref[0,]
# for(g in goi){
#   goi_ref = rbind(goi_ref,subset(ref, gene_name == g))
# }
ref_loci = GRanges(seqnames = ref$chrm, IRanges(ref$start, ref$end), gene_name = ref$gene_name, strand = ref$strand)
# ref_loci = subset(ref_loci, grepl("^HIST", gene_name))
# goi_loci = GRanges(seqnames = goi_ref$chrm, IRanges(goi_ref$start, goi_ref$end), gene_name = goi_ref$gene_name)
goi_loci = GRanges(seqnames = c("chr1", "chr6"), IRanges(c(149*10^6, 26*10^6), c(156*10^6, 28*10^6)))
# goi_loci = GRanges(seqnames = c("chr1"), IRanges(c(141*10^6), c(156*10^6)))

files = dir("bins_balanced/", full.names = T, pattern = "MCF10")
files= files[c(3,1)]

all_binned = list()
for(f in files){
  key = basename(f) %>% sub(pattern = ".save", replacement = "")
  print(f)
  load(f)
  all_binned[[key]] = as(binned, Class = "dgTMatrix")
}
names(all_binned) = c("MCF10A", "MCF10A-AT1")

my_plots = function(binned_name, loci_gr, gene_ref, ext){
  
  layout(1:3, heights = c(2,2,1.5))
  chrm = as.character(seqnames(loci_gr))
  
  start = start(loci_gr)
  start = as.integer((start - ext) / bin_size) * bin_size
  end = end(loci_gr)
  end = as.integer((end + ext + bin_size) / bin_size) * bin_size
  
  my_plotHic(binned_name = binned_name, 
             chrm = chrm, start = max(start,1), end = end,
             # loci_gr = loci_gr, 
             gene_ref = gene_ref)
  my_plotHic_lines(binned_name = binned_name,
                   chrm = chrm, start = max(start,1), end = end,
                   gene_ref = gene_ref,
                   min_dist = 4,
                   n_displayed = 100)
  #do bedgraph plot
  my_plotBdg(binned_name = binned_name, chrm = chrm, start = start, end = end)
  
}

if(exists("wd")){
  setwd(wd)
}else{
  wd = getwd()  
}

outdir = "TM_quick_checks"
dir.create(outdir, showWarnings = F)
outdir = normalizePath(outdir)
setwd(outdir)

goi_loci = GRanges("chr1", IRanges(173724478, 173924377))
goi_loci = GRanges("chr1", IRanges(173816175, 173877795))
for(i in (1:length(all_binned))){
  for(j in 1:length(goi_loci)){
    root = paste(sep = "_",
                 as.character(seqnames(goi_loci)[j]),
                 paste0(as.integer(mean(c(start(goi_loci)[j], end(goi_loci)[j]))/10^3), "K"),
                 names(all_binned)[i])
    print(root)
    # png(paste0(root, ".png"), width = 600, height = 800)
    pdf(paste0(root, ".pdf"), width = 6, height = 8)
    my_plots(binned_name = names(all_binned)[i], 
             loci_gr = goi_loci[j], 
             gene_ref = ref_loci, 
             ext = 40000)
    dev.off()
  }
}
setwd('..')


