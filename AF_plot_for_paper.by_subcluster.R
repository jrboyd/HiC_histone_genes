#uses HiC-Pro output
library(Sushi)
library(GenomicRanges)
library(magrittr)
library(data.table)
source("parse_gtf.R")
source("functions_my_HiC_plots.R")
library(ggplot2)

juice2triangle_plot = function(juice_file, bin_size = 40000){
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
j_files = dir("juicer_output/40kb/", pattern = "obs_vs_exp", full.names = T)

juice2binned = function(juice_file){
  dat = fread(juice_file)
  # dim_names = paste(chr, 1:(max(dat$V2, dat$V1) / bin_size))
  M = coord2Matrix(dat, bin_size = 40000, chr = "chr6")
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
full_ref_loci = GRanges(seqnames = ref$chrm, IRanges(ref$start, ref$end), gene_name = ref$gene_name, strand = ref$strand)

if(exists("wd")){
  setwd(wd)
}else{
  wd = getwd()  
}

outdir = "by_subcluster"
dir.create(outdir, showWarnings = F)
outdir = normalizePath(outdir)
setwd(outdir)

pdf("subclusters.wTRIM.pdf")
for(cl in names(all_binned)){
  binned = all_binned[[cl]]
  bin_size = 40000 
  bin_starts = sapply(strsplit(rownames(binned), " "), function(x){
    as.numeric(x[2]) * bin_size
  })
  
  
  m = as.matrix(binned); 
  # m = ifelse(m < 1, -1/m, m)
  m[is.infinite(m)] = 0
  
  #a large part of the matrix is 0, remove it
  k = range(which(colSums(m)> 0))
  k = k[1]:k[2]
  # heatmap.2(m[k, k], , trace = "n", Colv = F, Rowv = F, col = cr)
  
  bin_gr = GRanges(seqnames = "chr6", IRanges(bin_starts + 1, bin_starts + bin_size))
  
  my_hist_loci = full_ref_loci[grepl("HIST", full_ref_loci$gene_name)]
  hist_bin_hits = unique(queryHits(findOverlaps(bin_gr, my_hist_loci)))
  
  s = min(hist_bin_hits) - 25
  e = max(hist_bin_hits) + 40
  
  annotate_group = function(grp_regex, grp_color){
    grp_loci = full_ref_loci[grepl(grp_regex, full_ref_loci$gene_name)]
    grp_hits = queryHits(findOverlaps(bin_gr, grp_loci))
    grp_hits = intersect(grp_hits, names(grp_colors))
    grp_colors[as.character(grp_hits)] <<- grp_color
    grp_key_colors <<- c(grp_color, grp_key_colors)
    grp_names <<- c(gsub("[()^$]", "", grp_regex), grp_names)
  }
  
  grp_colors = rep("gray", length(s:e))
  names(grp_colors) = s:e
  grp_key_colors = "gray"
  grp_names = "other"
  ### add grp annotations here
  
  annotate_group("^HIST1", "green")
  # annotate_group("^HFE$", "orange")
  # annotate_group("^BTN", "purple")
  # annotate_group("^SLC", "darkgreen")
  # 
  # annotate_group("^OR", "red")
  # annotate_group("^ZNF", "darkblue")
  # annotate_group("^(Z|ZK)SCAN", "black")
  
  annotate_group("^TRIM", "blue")

  ###
  cs = which(!(grp_colors[-1] == grp_colors[-length(grp_colors)]))
  cr = rgb(colorRamp(c("white", "red"))(0:20/20)/255)
  #cap matrix values
  plot_m = m[s:e, s:e]
  raw_plot_m = plot_m
  m_max = quantile(plot_m, .999)
  plot_m = ifelse(plot_m > m_max, m_max, plot_m)
  plot_m = ifelse(plot_m < -m_max, -m_max, plot_m)
  #plot
  labs = sapply(strsplit(colnames(plot_m), " "), function(x){
    i = as.numeric(x[2])
    chr = x[1]
    start = i * 40000 + 1
    end = (i + 1) * 40000
    paste0(chr, ":", start,"-", end)
  })
  heatmap.2(plot_m, trace = "n", labRow = labs, labCol = labs, cexRow = .2, cexCol = .2,
            Colv = F, Rowv = F, 
            col = cr, 
            scale = "n",
            margins = c(8,8), 
            dendrogram = "n", 
            RowSideColors = grp_colors, ColSideColors = grp_colors, 
            rowsep = cs, colsep = cs, 
            sepcolor = "black", 
            main = paste(cl, "subclusters"), 
            key.title = "", key.xlab = "obs / exp", key.ylab = "", density.info = "n")
  legend("top", legend = grp_names, fill = grp_key_colors, bty = "n", ncol = 3)
  
  # heatmap.2(plot_m, trace = "n", Colv = F, Rowv = F, col = cr, margins = c(8,8), dendrogram = "n", RowSideColors = grp_colors, ColSideColors = grp_colors, rowsep = 1, colsep = cs, sepcolor = "black")
  
  starts = c(1, cs+1)
  ends = c(cs, nrow(plot_m))
  dn = paste(starts, ends, sep = "-")
  len = length(starts)
  sc_mat = matrix(0, len, len)
  sc_colors = character() 
  rownames(sc_mat) = dn
  colnames(sc_mat) = dn
  floor_m = raw_plot_m
  floor_m = ifelse(floor_m < 0, -1/floor_m, floor_m)
  for(i in 1:len){
    i_s = starts[i]
    i_e = ends[i]
    sc_colors = c(sc_colors, unique(grp_colors[i_s:i_e]))
    for(j in 1:len){
      j_s = starts[j]
      j_e = ends[j]
      block_m = floor_m[i_s:i_e, j_s:j_e, drop = F]
      #gets upper triangle of block
      # block_tri = sapply(1:ncol(block_m)-1, function(x){
      #   if(x > 0){
      #     diag(block_m[,-1:-x, drop = F])  
      #   }else{
      #     diag(block_m)
      #   }
      # })
      val = mean(block_m)
      sc_mat[i,j] = val
    }
  }
  cr2 = rgb(colorRamp(c("white", "red"))(0:20/20)/255)
  labs2 = apply(cbind(colnames(plot_m)[starts], colnames(plot_m)[ends]), 1, function(x){
    ssp = strsplit(x, " ")
    chr = ssp[[1]][1]
    start = as.numeric(ssp[[1]][2]) * 40000 + 1
    end = (as.numeric(ssp[[2]][2]) + 1) * 40000
    paste0(chr, ":", start,"-", end)
  })
  heatmap.2(sc_mat, col = cr2, labRow = labs2, labCol = labs, cexRow = 1.2, cexCol = 1.2,
            RowSideColors = sc_colors, 
            ColSideColors = sc_colors, 
            scale = "n",
            Rowv = F, Colv = F, 
            trace = "n", dendrogram = "n", 
            margins = c(14,14),
            main = paste(cl, "block averaged"), 
            key.title = "", key.xlab = "mean(obs / exp)", 
            key.ylab = "", density.info = "n")
  legend("top", legend = grp_names, fill = grp_key_colors, bty = "n", ncol = 3)
}
dev.off()
