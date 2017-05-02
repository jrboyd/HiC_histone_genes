#summarizes interchromosomal interactions focused on a particular region
library(magrittr)
library(GenomicRanges)
source("functions_my_HiC_plots.R")

load("hg38ref.save")
ref_gr = GRanges(seqnames = ref$chrm, IRanges(ref$start, ref$end), gene_name = ref$gene_name, gene_id = ref$gene_id)
goi = c("ZEB1", "PLAUR", "IL6", "CDH1", "FN1")
# goi = "FN1"
goi_ref = ref[0,]
for(g in goi){
  goi_ref = rbind(goi_ref,subset(ref, gene_name == g))
}

goi_loci = GRanges(seqnames = goi_ref$chrm, IRanges(goi_ref$start, goi_ref$end), gene_name = goi_ref$gene_name)


files = dir(path = "output_AF_hg38_integrated_filtered3_interchromosomal_v2/", pattern = "maskFull.save", full.names = T)
for(f in files){
  load(f)
  n_res = length(all_best)
  options(scipen = 90000000)
  all_best = lapply(all_best, function(x){
    chrms = sapply(strsplit(names(x), " "), function(x)x[1])
    i = as.numeric(sapply(strsplit(names(x), " "), function(x)x[2]))/40000
    best = numeric()
    for(chr in unique(chrms)){
      k = chrms == chr
      vals = x[k]
      pos = i[k]
      o = order(pos)
      vals = vals[o]
      pos = pos[o]
      adjacency = pos[-1] - pos[-length(pos)]
      kept = which(adjacency == 1)
      kept = sort(union(kept, kept+1))
      if(length(kept) == 0) next
      vals = vals[kept]
      names(vals) = paste(chr, pos[kept]*40000)
      best = c(best, vals)
    }
    return(best)
    # n = paste((sapply(strsplit(names(x), " "), function(x)x[1])),
    #       as.numeric(sapply(strsplit(names(x), " "), function(x)x[2])))
    # names(x) = n
    # return(x)
  })
  all_hits = unique(unlist(lapply(all_best, names)))
  m = matrix(0, nrow = length(all_hits), ncol = length(all_best))
  rownames(m) = all_hits
  colnames(m) = names(all_best)
  for(n in names(all_best)){
    best = all_best[[n]]
    m[names(best),n] = best
  }
  heatmap(m)
  
  
  
  df = as.data.frame(m)
  df$chrm = factor(sub("chr", "", sapply(strsplit(rownames(m), " "), function(x)x[1])))
  df$pos = as.integer(sapply(strsplit(rownames(m), " "), function(x)x[2]))
  df$i = df$pos / 40000 
  
  o = order(df$i)
  df = df[o,]
  o = order(as.integer(as.character(df$chrm)))
  df = df[o,]
  
  # df[,1:5] = apply(df[,1:5], 2, function(x)x/max(x))
  # heatmap.2(as.matrix(df[,1:5]), Colv = T, Rowv = F, trace = "n")
  
  
  
  head(df)
  all_gene_interaction = list()
  for(i in 1:n_res){
    o = order(df[,i], decreasing = T)
    rn = rownames(df)[o]
    chrm = sapply(strsplit(rn, " "), function(x)x[1])
    start = as.integer(sapply(strsplit(rn, " "), function(x)x[2]))
    end = start + 40000
    hic_gr = GRanges(chrm, IRanges(start, end), values = df[rn,i])
    olaps = findOverlaps(query = hic_gr, subject = ref_gr)
    gene_interactions = data.frame(gene_name = ref_gr[subjectHits(olaps)]$gene_name, microinteraction = hic_gr[queryHits(olaps)]$values * 10^6)
    all_gene_interaction[[colnames(df)[i]]] = gene_interactions
  }
  all_g = unique(as.character(unlist(lapply(all_gene_interaction, function(x)x[,1]))))
  gene_interactions_mat = matrix(0, nrow = length(all_g), ncol = n_res)
  colnames(gene_interactions_mat) = names(all_gene_interaction)
  rownames(gene_interactions_mat) = all_g
  for(i in 1:length(all_gene_interaction)){
    agi = all_gene_interaction[[i]]
    gene_interactions_mat[as.character(agi[,1]),i] = agi[,2]
  }
  
}

