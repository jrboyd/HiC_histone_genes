#uses HiC-Pro output
library(Sushi)
library(GenomicRanges)
library(magrittr)
library(data.table)
source("parse_gtf.R")
source("functions_my_HiC_plots.R")
library(ggplot2)

juice2triangle_plot = function(juice_file, bin_size = 50000){
  # juice_file = "~/juicebox/mcf10a_matrix_obs_vs_exp.matrix"
  dat = read.table(juice_file)
  colnames(dat) = c("i_gen_start", "j_gen_start", "val")
  dat$i_gen_end = dat$i_gen_start + bin_size
  dat$j_gen_end = dat$j_gen_start + bin_size
  start = 25700000
  end = 29100000 - bin_size
  rotate_xy = function(x, y, theta){
    rot = theta/360  * 2* pi
    rot_mat = matrix(c(cos(rot), sin(rot), -sin(rot), cos(rot)), ncol = 2)
    my_points = cbind(x, y)
    rot_points = t(crossprod(rot_mat, t(my_points)))
    rot_points
  }
  xy = rotate_xy(dat$i_gen_start)
  
  ggplot(subset(dat, i_gen_start >= start & i_gen_start <= end & j_gen_start >= start & j_gen_start <= end
                )) + geom_tile(aes(x = i_gen_start, y = j_gen_start,  fill = val)) + coord_fixed() 
}


coord2Matrix = function(x, bin_size = 50000, chr = "chr6"){
  x[, V1 := V1 / bin_size]
  x[, V2 := V2 / bin_size]
  # x = x[V1 %in% to_keep & V2 %in% to_keep]
  # x[, c("V1", "V2") := .(V1 - min(to_keep) , V2 - min(to_keep) )]
  diag = x[V1 == V2]
  main = x[V1 != V2]
  x = rbind(diag, main, main[, .(V1 = V2, V2 = V1, V3 = V3)])
  len = max(x$V1, x$V2) + 5
  M = new("dgTMatrix",
          i = as.integer(x$V1-1),
          j = as.integer(x$V2-1), 
          x=x$V3, 
          Dim= as.integer(c(len, len)))
  # m = as.matrix(M)
  rownames(M) = paste(chr, 1:nrow(M))
  colnames(M) = paste(chr, 1:ncol(M))
  M
}

# j_files = dir("~/juicebox/", pattern = "obs2", full.names = T)
j_files = dir("juicer_output/40kb", pattern = "obs_vs_exp", full.names = T)
# j_files = dir("juicer_output/40kb", pattern = "obs_40", full.names = T)
bin_size = 40000

juice2binned = function(juice_file){
  dat = fread(juice_file)
  # dim_names = paste(chr, 1:(max(dat$V2, dat$V1) / bin_size))
  M = coord2Matrix(dat, bin_size = bin_size, chr = "chr6")
  M
}

all_binned = lapply(j_files, juice2binned)
names(all_binned) = c("MCF10A", "MCF10A-AT1", "MCF10A-CA1a")

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

#reduce sytenic regions of similarly prefixed gene names to 1 entry per continuous group
simplify_ref = function(in_ref, target = "^HIST1"){
  #matches for target
  hits = which(grepl(target, in_ref$gene_name))
  out_ref = in_ref[-hits,]
  ranges = matrix(ncol = 2, nrow = 0)
  while(length(hits) > 0){
    s = hits[1]
    hits = hits[-1]
    e = s
    for(i in hits){
      if(i == e + 1){
        e = i
        hits = hits[-1]
      }else{
        break
      }
    }
    ranges = rbind(ranges, cbind(s, e))
  }
  for(i in 1:nrow(ranges)){
    x = ranges[i,]
    n = x[2] - x[1] + 1
    chrm = unique(in_ref[x,]$chrm)
    start = min(in_ref[x,]$start)
    end = max(in_ref[x,]$end)
    name = paste(sub("[^$]", "", target), "c", i, "s", n)
    to_add = data.frame(gene_id = name, gene_name = name, chrm = chrm, start = start, end = end, strand = "*", gene_type = "cluster" )
    out_ref  = rbind(out_ref, to_add)
  }
  out_ref = out_ref[order(out_ref$start),]
  out_ref = out_ref[order(out_ref$chrm),]
  out_ref
}
simp_ref = ref
simp_ref = simplify_ref(simp_ref, "^HIST1")
simp_ref = simplify_ref(simp_ref, "^ZNF")
simp_ref = simplify_ref(simp_ref, "^BTN")
simp_ref = simplify_ref(simp_ref, "^(Z|ZK)SCAN")
simp_ref = simplify_ref(simp_ref, "^SLC")
simp_ref = simplify_ref(simp_ref, "^OR")
simp_ref = simplify_ref(simp_ref, "^LINC")
ref = simp_ref
ref_loci = GRanges(seqnames = ref$chrm, IRanges(ref$start, ref$end), gene_name = ref$gene_name, strand = ref$strand)
goi_loci = GRanges(seqnames = c("chr1", "chr6"), IRanges(c(149*10^6, 24750000), c(156*10^6, 30250000)))
my_loci = goi_loci[2]

my_plots = function(binned_name, loci_gr, gene_ref, ext, quantile_cap = 1, zrange = NULL){
  
  layout(1:3, heights = c(2,2,1.5))
  chrm = as.character(seqnames(loci_gr))
  
  start = start(loci_gr)
  start = as.integer((start - ext) / bin_size) * bin_size
  end = end(loci_gr)
  end = as.integer((end + ext + bin_size) / bin_size) * bin_size
  print("triangle plot")
  my_plotHic(binned_name = binned_name, 
             chrm = chrm, start = max(start,1), end = end, applyLog = F,
             # loci_gr = loci_gr, 
             gene_ref = gene_ref, quantile_cap = quantile_cap, zrange = zrange)
  print("lines plots")
  my_plotHic_lines(binned_name = binned_name,
                   chrm = chrm, start = max(start,1), end = end,
                   gene_ref = gene_ref, bin_size = bin_size,
                   min_dist = 4,
                   n_displayed = 100, quantile_cap = 1, zrange = zrange)
  #do bedgraph plot
  print("bdg plots")
  my_plotBdg(binned_name = binned_name, chrm = chrm, start = start, end = end)
  
}

if(exists("wd")){
  setwd(wd)
}else{
  wd = getwd()  
}

outdir = "output_AF_paper_obs_vs_exp_40kb"
dir.create(outdir, showWarnings = F)
outdir = normalizePath(outdir)
setwd(outdir)


for(i in (1:length(all_binned))){
  for(j in 1:length(my_loci)){
    root = paste(sep = "_",
                 as.character(seqnames(my_loci)[j]),
                 paste0(as.integer(mean(c(start(my_loci)[j], end(my_loci)[j]))/10^3), "K"),
                 names(all_binned)[i])
    print(root)
    # png(paste0(root, ".png"), width = 600, height = 800)
    pdf(paste0(root, "_obsvsexp_40kb_scale_quantile.pdf"), width = 6, height = 8)
    my_plots(binned_name = names(all_binned)[i],
             loci_gr = my_loci[j], 
             gene_ref = ref_loci, 
             ext = 0, quantile_cap = .999, zrange = NULL)
    dev.off()
  }
}
setwd('..')




