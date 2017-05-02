#summarizes interchromosomal interactions focused on a particular region
library(magrittr)
library(GenomicRanges)
source("functions_my_HiC_plots.R")

href_file = "hist_ref.save"
if(file.exists(href_file)){
  load(href_file)
}else{
  source('parse_gtf.R')
  ref = parse_gtf("gencode.v24.annotation.gtf", rownames_attrib = "gene_id", feature_type = "gene")
  is_hist = grepl("^HIST", ref$gene_name)
  hist_ref = ref[is_hist,]
  ref[ref$gene_name == "NPAT",]
  hist_ref = rbind(hist_ref, ref[ref$gene_name == "NPAT",], ref[ref$gene_name == "HINFP" ,])
  save(ref, hist_ref, file = href_file)
}

ref = ref[!duplicated(ref$gene_name),]
ref_gr = GRanges(ref$chrm, IRanges(ref$start, ref$end), gene_id = ref$gene_id, gene_name = ref$gene_name)

loci_gap = 10^6
hist_pos = GRanges(seqnames = hist_ref$chrm, IRanges(hist_ref$start - loci_gap, hist_ref$end + loci_gap))
names(hist_pos) = hist_ref$gene_id
hist_loci = reduce(hist_pos)



files = dir(path = "interChr_results/", pattern = "Intra.save", full.names = T)
files = dir(path = "interChr_results/", pattern = "2MB.save", full.names = T)
for(f in files[10]){
  load(f)
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
  for(i in 1:5){
    o = order(df[,i], decreasing = T)
    rn = rownames(df)[o][1:5]
    chrm = sapply(strsplit(rn, " "), function(x)x[1])
    start = as.integer(sapply(strsplit(rn, " "), function(x)x[2]))
    end = start + 40000
    hic_gr = GRanges(chrm, IRanges(start, end), values = df[rn,i])
    olaps = findOverlaps(query = hic_gr, subject = ref_gr)
    gene_interactions = data.frame(gene_name = ref_gr[subjectHits(olaps)]$gene_name, microinteraction = hic_gr[queryHits(olaps)]$values * 10^6)
    all_gene_interaction[[colnames(df)[i]]] = gene_interactions
  }
  all_g = unique(as.character(unlist(lapply(all_gene_interaction, function(x)x[,1]))))
  gene_interactions_mat = matrix(0, nrow = length(all_g), ncol = 5)
  colnames(gene_interactions_mat) = names(all_gene_interaction)
  rownames(gene_interactions_mat) = all_g
  for(i in 1:length(all_gene_interaction)){
    agi = all_gene_interaction[[i]]
    gene_interactions_mat[as.character(agi[,1]),i] = agi[,2]
  }
  
}

