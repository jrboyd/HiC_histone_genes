library(Sushi)
library(GenomicRanges)
library(magrittr)
source("H:/R_workspace/jrb_R_scripts/parse_gtf.R")
source("functions_my_HiC_plots.R")
href_file = "hist_ref_19.save"
if(file.exists(href_file)){
  load(href_file)
}else{
  ref = parse_gtf("C:/Users/jrboyd/Downloads/gencode.v19.annotation.gtf", rownames_attrib = "gene_id", feature_type = "gene")
  is_hist = grepl("^HIST", ref$gene_name)
  hist_ref = ref[is_hist,]
  ref[ref$gene_name == "NPAT",]
  hist_ref = rbind(hist_ref, ref[ref$gene_name == "NPAT",], ref[ref$gene_name == "HINFP" ,])
  save(hist_ref, file = href_file)
}
files = c(
  "balanced_bins/Galaxy4-[hESC_TotUnbalMatrix_hg19_40kb].data",
  "balanced_bins/Galaxy2-[10A_TotUnbalMatrix_hg19_40kb].data",
  "balanced_bins/Galaxy3-[AT1_TotUnbalMatrix_hg19_40kb].data",
  "balanced_bins/Galaxy1-[CA1a_TotUnbalMatrix_hg19_40kb].data",
  "balanced_bins/Galaxy5-[MCF7_TotUnbalMatrix_hg19_40kb].data"
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

loci_gap = 10^6
hist_pos = GRanges(seqnames = hist_ref$chrm, IRanges(hist_ref$start - loci_gap, hist_ref$end + loci_gap))
names(hist_pos) = hist_ref$gene_id
hist_loci = reduce(hist_pos)

outdir = "output_hg19unbalanced_hmaps"
dir.create(outdir)
setwd(outdir)
ext = 3*10^6
for(i in 1:length(all_binned)){
  for(j in 1:length(hist_loci)){
    root = paste(sep = "_",
                 as.character(seqnames(hist_loci)[j]),
                 paste0(as.integer(mean(c(start(hist_loci)[j], end(hist_loci)[j]))/10^3), "K"),
                 names(all_binned)[i])
    print(root)
    png(paste0(root, ".png"), width = 800, height = 600)
    
    my_plotHic(names(all_binned)[i], loci_gr = hist_loci[j], hist_ref, p_ext = ext)
    dev.off()
  }
  print("NPAT and HINFP")
  loci = hist_loci[length(hist_loci)-1]
  end(loci) = end(hist_loci[length(hist_loci)])
  png(paste0("NPAT-and-HINFP_", names(all_binned)[i], ".png"), width = 800, height = 600)
  my_plotHic(binned_name = names(all_binned)[i], loci_gr = loci, gene_ref = hist_ref, p_ext = ext)
  dev.off()
}
setwd('..')

outdir = "output_hg19unbalanced_lines"
dir.create(outdir)
setwd(outdir)
for(i in 1:length(all_binned)){
  for(j in 1:length(hist_loci)){
    root = paste(sep = "_",
                 as.character(seqnames(hist_loci)[j]),
                 paste0(as.integer(mean(c(start(hist_loci)[j], end(hist_loci)[j]))/10^3), "K"),
                 names(all_binned)[i])
    print(root)
    png(paste0(root, ".png"), width = 800, height = 600)
    chrm = as.character(seqnames(hist_loci)[j])
    start = start(hist_loci)[j]
    end = end(hist_loci)[j]
    my_plotHic_lines(binned_name = names(all_binned)[i],
                     chrm = chrm, start = start - ext, end = end + ext,
                     gene_ref = hist_ref,
                     min_dist = 4,
                     n_displayed = 200)
    dev.off()
  }
  print("NPAT and HINFP")
  loci = hist_loci[length(hist_loci)-1]
  end(loci) = end(hist_loci[length(hist_loci)])
  png(paste0("NPAT-and-HINFP_", names(all_binned)[i], ".png"), width = 800, height = 600)
  chrm = as.character(seqnames(loci))
  start = start(loci)
  end = end(loci)
  my_plotHic_lines(binned_name = names(all_binned)[i],
                   chrm = chrm, start = start - ext, end = end + ext,
                   gene_ref = hist_ref,
                   min_dist = 10,
                   n_displayed = 200)
  dev.off()
}
setwd('..')
