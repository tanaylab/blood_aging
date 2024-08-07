library(pheatmap)
library(zoo)
library(dplyr)
library(metacell)
library(openxlsx)
library(tgstat)
library(Matrix)

get.one.samp.per.indiv <- function(all.samp, our.cdata, cell.to.samp) {
  samp.to.indiv = sapply(strsplit(all.samp, '-'), function(x) x[2])
  names(samp.to.indiv) = all.samp

  stopifnot(length(cell.to.samp) == nrow(our.cdata$obs))
  stopifnot(all(names(cell.to.samp) == rownames(our.cdata$obs)))

  samp.src.tbl = table(cell.to.samp, as.character(our.cdata$obs$src))
  stopifnot(all(rowSums(samp.src.tbl > 0) == 1))
  samp.to.src = colnames(samp.src.tbl)[apply(samp.src.tbl, 1, which.max)]
  names(samp.to.src) = rownames(samp.src.tbl)

  samp.tech.tbl = table(cell.to.samp, as.character(our.cdata$obs$is_ultima))
  stopifnot(all(rowSums(samp.tech.tbl > 0) == 1))
  samp.to.is.ultima = as.logical(colnames(samp.tech.tbl)[apply(samp.tech.tbl, 1, which.max)])
  names(samp.to.is.ultima) = rownames(samp.tech.tbl)

  all.indivs = unique(samp.to.indiv)
  indiv.to.sel.samp = sapply(all.indivs, function(cur.indiv) {
    all.indiv.samp = names(samp.to.indiv)[samp.to.indiv == cur.indiv]
    if (length(all.indiv.samp) == 1) {
      return(all.indiv.samp[1])
    }

    # explicit mapping for some individuals for which the ultima library was better than the illumina library
    # because it was sequenced more deeply, and the illumina cut away some cells, so less cells were cut in ultima
    explicit.samps = list("N152"="demux_22_02_21_ultima-N152", "N153"="demux_22_02_21_ultima-N153", "N155"="demux_22_02_21_ultima-N155", 
                          "N156"="demux_22_02_21_ultima-N156", "N159"="demux_22_02_21_ultima-N159",
	                  "N91"="demux_01_02_21_ultima-N91", "N130"="demux_01_02_21_ultima-N130", 
	                  "N133"="demux_01_02_21_ultima-N133", "N135"="demux_01_02_21_ultima-N135", 
	                  "N166"="demux_07_03_21_ultima-N166", "N167"="demux_07_03_21_ultima-N167", "N168"="demux_07_03_21_ultima-N168", 
	                  "N169"="demux_07_03_21_ultima-N169", "N171"="demux_07_03_21_ultima-N171", "N176"="demux_07_03_21_ultima-N176")

    samp.srcs = samp.to.src[all.indiv.samp]
    samp.is.ultima = samp.to.is.ultima[all.indiv.samp]
    if (cur.indiv %in% names(explicit.samps)) {
      stopifnot(explicit.samps[[cur.indiv]] %in% all.indiv.samp)
      return(explicit.samps[[cur.indiv]])
    }
    if (sum(samp.srcs != 'new') == 1) {
      return(all.indiv.samp[which(samp.srcs != 'new')])
    }
    if (sum(!samp.is.ultima) == 1) {
      return(all.indiv.samp[which(!samp.is.ultima)])
    }
    indiv.to.used.samp.map = list("N180" = "demux_10_01_22-N180", 
                                  "N307" = "demux_04_12_22-N307")

    stopifnot(cur.indiv %in% names(indiv.to.used.samp.map))
    stopifnot(indiv.to.used.samp.map[cur.indiv] %in% all.indiv.samp)
    return(indiv.to.used.samp.map[[cur.indiv]])
  })
  return(indiv.to.sel.samp)

}

plot.age.distributions <- function(our.cdata, our.annotation, our.assign, fig.dir) {
  cell.to.meas = sprintf('%s-%s', our.cdata$obs$exp_name, our.cdata$obs$indiv_id)
  names(cell.to.meas) = rownames(our.cdata$obs)
  cell.to.samp = gsub('_\\d_ultima-', '_ultima-', gsub('_\\d-', '-', cell.to.meas))
  tmp.tbl = table(cell.to.samp, our.annotation[our.assign[names(cell.to.samp)]])
  all.samp = unique(cell.to.samp)

  indiv.to.used.samp = get.one.samp.per.indiv(all.samp, our.cdata, cell.to.samp)

  age.and.sex = get.age.and.sex(our.cdata, indiv.to.used.samp, cell.to.samp)
  ages = age.and.sex$ages
  sexes = age.and.sex$sexes

  all.indivs = as.character(unique(our.cdata$obs$indiv_id))
  is.male = sexes == 'male'
  all.males = all.indivs[sexes[all.indivs] == 'male']
  all.females = all.indivs[sexes[all.indivs] == 'female']
  indivs.cols = ifelse(is.male[all.indivs], 'lightgreen', '#FFD580')

  png(file.path(fig.dir, 'age_distributions.png'), w=800, h=400)
  par(mfrow=c(1, 2))
  hist(ages[all.males], col=most.common(indivs.cols[all.males]), xlim=range(ages[all.indivs]) + c(-5, 10))
  hist(ages[all.females], col=most.common(indivs.cols[all.females]), xlim=range(ages[all.indivs]) + c(-5, 10))
  dev.off()
}

fig1 <- function() {
  dir.create(BASE.FIG.DIR, showWarnings=F)
  fig.dir = file.path(BASE.FIG.DIR, 'fig1')
  dir.create(fig.dir, showWarnings=F)

  our.cdata = anndata::read_h5ad(file.path(MODEL.DIR, '148_indiv_ref_cells.h5ad'))
  our.mdata = anndata::read_h5ad(file.path(MODEL.DIR, '148_indiv_ref_metacells.h5ad'))
  tmp.umi.counts = our.mdata$layers[['total_umis']]
  our.legc = log2(1e-5 + t(tmp.umi.counts / rowSums(tmp.umi.counts)))
  colnames(our.legc) = paste0('mc', 0:(ncol(our.legc) - 1))
  our.assign = paste0('mc', our.cdata$obs$metacell)
  names(our.assign) = rownames(our.cdata$obs)

  our.annotation = annotate.model(our.legc)

  # distribution across metacells
  cur.used.ctypes = c('purple', 'red', 'plum', '#EEBB6E', "#7F9D00", 'gold', 'brown', "#6E1C71", 'blue', 'darkblue', 'lightblue')
  used.mcs = names(our.annotation)[our.annotation %in% cur.used.ctypes]
  used.cells = names(our.assign)[our.assign %in% used.mcs]
  indiv.mc.tbl = table(as.character(our.cdata$obs[used.cells, 'indiv_id']), our.assign[used.cells])
  print(quantile(colSums(indiv.mc.tbl > 0)))

  plot.age.distributions(our.cdata, our.annotation, our.assign, fig.dir)

  full.umap = plot.full.umap(our.legc, our.annotation, fig.dir)
  filtered.umap = plot.filtered.umap(our.legc, our.annotation, fig.dir)

  # main gene scatters
  good.mcs = names(our.annotation)[our.annotation != 'grey']
  png(file.path(fig.dir, 'cd34_vs_avp.png'), w=400, h=400)
  plot(our.legc['CD34', good.mcs], our.legc['AVP', good.mcs], pch=21, bg=our.annotation[good.mcs], xlab='CD34', ylab='AVP')
  abline(v=-14.5,lwd=2)
  dev.off()

  png(file.path(fig.dir, 'clp_scatter.png'), width=600, height=600)
  plot(our.legc['RUNX3', used.mcs], our.legc['DNTT', used.mcs], pch=21, bg=our.annotation[used.mcs], cex=1.5)
  dev.off()

  png(file.path(fig.dir, 'gmp_scatter.png'), width=600, height=600)
  plot(our.legc['GATA1', used.mcs] - our.legc['VPREB1', used.mcs], our.legc['MPO', used.mcs], pch=21, bg=our.annotation[used.mcs], cex=1.5)
  dev.off()

  # gene x cell type expression heatmap
  ctype.to.col = c("CLP-E"="#6E1C71", "GMP-E"="#7F9D00", "MEBEMP-M"="#EEBB6E", "Endothel"="black", "CLP-M"="blue", "HSC"="brown", 
                   "Monocytes"="chartreuse", "DC"="cyan", "CLP-L"="darkblue", "GMP-L"="darkgreen", "MPP"="gold", "NKTDP"="lightblue", "NKT"="orange", 
		   "MEBEMP-L"="plum", "BEMP"="purple", "EP"="red", "B"="steelblue", "Unknown"="white")
  col.to.ctype = names(ctype.to.col)
  names(col.to.ctype) = ctype.to.col
  selected.genes = c('LMO4', 'MS4A2', 'HBD', 'APOC1', 'KLF1', 'GATA1', 'GATA2', 'AVP', 'HLF', 'MPO', 'AZU1', 
                     'HOPX', 'VPREB1', 'DNTT', 'RUNX3', 'ID2', 'ACY3', 'IRF8', 'SPIB', 'CD74', 'CD79A', 'CD79B', 
		     'S100A9', 'LYZ', 'CD3D', 'GNLY', 'IL7R', 'CAV1', 'ADIRF', 'CD34')
  #color.order = c('purple', 'red', 'plum', '#EEBB6E', 'gold', 'brown', '#7F9D00', 'darkgreen', '#6E1C71', 
  #                'blue', 'darkblue', '#566CF2', 'lightblue', 'cyan', 'steelblue', 'chartreuse', 'orange', 'black')
  color.order = c('purple', 'red', 'plum', '#EEBB6E', 'gold', 'darkgreen', '#7F9D00', 'brown', '#6E1C71', 
                  'blue', 'darkblue', 'lightblue', 'cyan', 'steelblue', 'chartreuse', 'orange', 'black')
  our.annotation.tmp = our.annotation[!(our.annotation %in% c('grey', 'white'))]
  mean.per.type = do.call(cbind, tapply(names(our.annotation.tmp), our.annotation.tmp, function(mc.names) rowMeans(our.legc[,mc.names, drop=F])))
  mean.per.type.norm = mean.per.type - apply(mean.per.type, 1, median)
  tmp.df = data.frame(mc_col=col.to.ctype[color.order])
  rownames(tmp.df) = col.to.ctype[color.order]
  tmp.colors = unique(our.annotation)
  names(tmp.colors) = tmp.colors
  names(tmp.colors) = col.to.ctype[tmp.colors]

  shades = colorRampPalette(c('darkblue', 'white', 'darkred'))(100)
  mean.per.type.norm.ord = mean.per.type.norm[selected.genes, color.order]
  colnames(mean.per.type.norm.ord) = col.to.ctype[colnames(mean.per.type.norm.ord)]
  pheatmap(pmin(pmax(mean.per.type.norm.ord, -4), 4), annotation_col=tmp.df, annotation_colors=list(mc_col=tmp.colors), col=shades, breaks=seq(-4, 4, length.out=100),
           cluster_rows=F, cluster_cols=F, show_colnames=F, filename=file.path(fig.dir, 'marker_heatmap.png'), height=10, width=8)


  plot.hsc.modes.and.dynamics(our.legc, our.annotation, fig.dir=fig.dir)

  mk.model.analysis(fig.dir=fig.dir)
  bm.comparison(our.mdata, our.annotation, fig.dir=fig.dir)
  bm.comparison.nktdp(our.legc, fig.dir=fig.dir)

  # scatters introducing the NKTDP and early BEMP populations
  png(file.path(fig.dir, "hbd_lmo4.png"), w=400,h=400)
  f_filt = !(our.annotation %in% c("orange", "cyan","chartreuse", "steelblue", "darkgreen", "black", "grey"))
  plot(our.legc['HBD', f_filt], our.legc['LMO4', f_filt], pch=19, col=our.annotation[f_filt], xlab='HBD', ylab='LMO4')
  dev.off()
  
  png(file.path(fig.dir, "irf8_tcf7.png"), w=400,h=400)
  f_filt = !our.annotation %in% setdiff(unique(our.annotation), c("darkblue","lightblue"))
  plot(our.legc['IRF8', f_filt], our.legc['TCF7', f_filt], pch=19, col=our.annotation[f_filt], xlab='IRF8', ylab='TCF7')
  dev.off()

  bemp.nktdp.scatters(our.legc, our.annotation, fig.dir)
}

bm.comparison <- function(our.mdata, our.annotation, fig.dir=file.path(BASE.FIG.DIR, 'fig1')) {
  bm.atlas.mdata = anndata::read_h5ad(BM.MODEL.PATH)
  bm.atlas.legc = t(log2(1e-5 + bm.atlas.mdata$X / rowSums(bm.atlas.mdata$X)))
  bm.annotation.type = as.character(bm.atlas.mdata$obs$type)
  bm.colors = read.csv(BM.MODEL.ANNOTATION, stringsAsFactors=F)
  rownames(bm.colors) = bm.colors[,1]
  bm.annotation = bm.colors[bm.annotation.type, 2]
  names(bm.annotation) = rownames(bm.atlas.mdata$obs)

  # load setty et al data (from palantir paper)
  pal.cdata = anndata::read_h5ad(file.path(MODEL.DIR, 'bm_palantir_ref_cells.h5ad'))
  pal.mdata = anndata::read_h5ad(file.path(MODEL.DIR, 'bm_palantir_ref_metacells.h5ad'))
  pal.legc = t(log2(1e-5 + pal.mdata$X / rowSums(pal.mdata$X)))

  our.bm.mdata = anndata::read_h5ad(file.path(MODEL.DIR, 'bm_our_ref_metacells.h5ad'))
  our.bm.mdata.x = our.bm.mdata$layers[['total_umis']]
  our.bm.legc = data.matrix(t(log2(1e-5 + our.bm.mdata.x / rowSums(our.bm.mdata.x))))
  colnames(our.bm.legc) = paste0('mc', 0:(ncol(our.bm.legc) - 1))

  s.genes = c('TYMS', 'H2AFZ', 'PCNA', 'MCM4', 'HELLS', 'MKI67')
  pb.ctypes.ord.for.s = c('purple', 'red', 'plum', '#EEBB6E', '#7F9D00', 'gold', 'brown', '#6E1C71', 'blue', 'darkblue', 'lightblue')
  bm.ctypes.ord.for.s = c('purple', 'red', 'plum', 'darkgreen', 'gold', 'azure4', 'blue', 'lightblue')
  pb.s.scores.by.ctype = split(colMeans(our.legc[s.genes, names(our.annotation)]), our.annotation)[pb.ctypes.ord.for.s]
  bm.s.scores.by.ctype = split(colMeans(bm.atlas.legc[s.genes, names(bm.annotation)]), bm.annotation)[bm.ctypes.ord.for.s]
  s.ylim = range(c(unlist(pb.s.scores.by.ctype), unlist(bm.s.scores.by.ctype)))
  png(file.path(fig.dir, 'bm_s_scores.png'))
  boxplot(bm.s.scores.by.ctype, col=names(bm.s.scores.by.ctype), ylim=s.ylim, las=2)
  dev.off()
  png(file.path(fig.dir, 'our_s_scores.png'))
  boxplot(pb.s.scores.by.ctype, col=names(pb.s.scores.by.ctype), ylim=s.ylim, las=2)
  dev.off()
  png(file.path(fig.dir, 'bm_and_pb_s_scores.png'))
  tmp.s.list = c(pb.s.scores.by.ctype, list(NA, NA, NA), bm.s.scores.by.ctype)
  boxplot(tmp.s.list, col=c(names(pb.s.scores.by.ctype), NA, NA, NA, names(bm.s.scores.by.ctype)), ylim=s.ylim, las=2)
  dev.off()

  # our model uses "lateral_gene", and palantir's "forbidden" because they were built with different mc versions
  our.bad.genes = rownames(our.mdata$var)[our.mdata$var$lateral_gene]
  pal.bad.genes = rownames(pal.mdata$var)[pal.mdata$var$forbidden_gene]
  genes.for.proj = intersect(setdiff(rownames(bm.atlas.mdata$var)[bm.atlas.mdata$var$top_feature_gene], c(our.bad.genes, pal.bad.genes)),
                             intersect(rownames(our.legc), rownames(pal.legc)))
  tmp.our.cors = cor(bm.atlas.legc[genes.for.proj,], our.legc[genes.for.proj,])
  our.top5.matches = lapply(1:ncol(tmp.our.cors), function(j) which(tmp.our.cors[,j] >= sort(tmp.our.cors[,j], decreasing=T)[5]))
  our.proj.annotation = sapply(our.top5.matches, function(x) most.common(bm.annotation[x]))
  names(our.proj.annotation) = colnames(our.legc)

  png(file.path(fig.dir, 'bm_umap.png'), width=600, height=600)
  plot(bm.atlas.mdata$obs$umap_x, bm.atlas.mdata$obs$umap_y, bg=bm.annotation[rownames(bm.atlas.mdata$obs)], pch=21, yaxt='n', xaxt='n')
  dev.off()
  bm.colors.ord = c('gold', 'azure4', 'plum', 'pink', 'red', 'purple', 
                    'darkgreen', 'darkseagreen2', 'chartreuse', 
		    'lightblue', 'cyan', 'blue', 'steelblue4', 'steelblue', 'cadetblue', 'orange')
  stopifnot(sort(unique(bm.annotation)) == sort(bm.colors.ord))
  png(file.path(fig.dir, 'bm_annotation_legend.png'), width=500, height=2500)
  plot(rep(1, length(bm.colors.ord)), (1:length(bm.colors.ord)) / 3, pch=19, col=bm.colors.ord, cex=10)
  dev.off()

  umap.coords = data.frame(bm.atlas.mdata$obs$umap_x, bm.atlas.mdata$obs$umap_y)
  rownames(umap.coords) = rownames(bm.atlas.mdata$obs)
  projected.coords = do.call(rbind, lapply(1:length(our.top5.matches), function(i) colMeans(umap.coords[colnames(bm.atlas.legc)[our.top5.matches[[i]]],])))
  rownames(projected.coords) = colnames(our.legc)
  #plot(bm.atlas.mdata$obs$umap_x, bm.atlas.mdata$obs$umap_y, bg=bm.annotation[rownames(bm.atlas.mdata$obs)], pch=21)
  #points(projected.coords[,1], projected.coords[,2], bg='white', pch=21)
  set.seed(42)
  umap.noise.x = rnorm(nrow(projected.coords), 0, 0.3)
  umap.noise.y = rnorm(nrow(projected.coords), 0, 0.3)
  names(umap.noise.x) = rownames(projected.coords)
  names(umap.noise.y) = rownames(projected.coords)
  # some B metacells have ambient high CD34
  mcs.to.project = setdiff(colnames(our.legc)[our.legc['CD34',] > -14.5], names(our.annotation)[our.annotation %in% c('grey', 'steelblue')])
  png(file.path(fig.dir, 'our_and_bm_atlas_2d.png'), height=500, width=500)
  plot(bm.atlas.mdata$obs$umap_x, bm.atlas.mdata$obs$umap_y, bg='white', pch=21, cex=0.8, yaxt='n', xaxt='n')
  points(projected.coords[mcs.to.project,1] + umap.noise.x[mcs.to.project], 
         projected.coords[mcs.to.project,2] + umap.noise.y[mcs.to.project], bg=our.annotation[mcs.to.project], pch=21, cex=1)
  dev.off()

  # project pal
  tmp.pal.cors = cor(bm.atlas.legc[genes.for.proj,], pal.legc[genes.for.proj,])
  pal.top5.matches = lapply(1:ncol(tmp.pal.cors), function(j) which(tmp.pal.cors[,j] >= sort(tmp.pal.cors[,j], decreasing=T)[5]))
  pal.proj.annotation = sapply(pal.top5.matches, function(x) most.common(bm.annotation[x]))
  names(pal.proj.annotation) = colnames(pal.legc)

  pal.projected.coords = do.call(rbind, lapply(1:length(pal.top5.matches), function(i) colMeans(umap.coords[colnames(bm.atlas.legc)[pal.top5.matches[[i]]],])))
  rownames(pal.projected.coords) = colnames(pal.legc)
  set.seed(43)
  umap.noise.x = rnorm(nrow(pal.projected.coords), 0, 0.3)
  umap.noise.y = rnorm(nrow(pal.projected.coords), 0, 0.3)
  names(umap.noise.x) = rownames(pal.projected.coords)
  names(umap.noise.y) = rownames(pal.projected.coords)
  pal.mcs.to.project = colnames(pal.legc)[pal.legc['CD34',] > -14.5] 
  png(file.path(fig.dir, 'pal_and_bm_atlas_2d.png'), height=500, width=500)
  plot(bm.atlas.mdata$obs$umap_x, bm.atlas.mdata$obs$umap_y, bg='white', pch=21, cex=0.8, yaxt='n', xaxt='n')
  points(pal.projected.coords[pal.mcs.to.project,1] + umap.noise.x[pal.mcs.to.project], 
         pal.projected.coords[pal.mcs.to.project,2] + umap.noise.y[pal.mcs.to.project], bg=pal.proj.annotation[pal.mcs.to.project], pch=21, cex=1)
  dev.off()

  # project our bm data
  our.bm.dbl.mcs = colnames(our.legc)[(our.bm.legc['LMO4',] > -11.5 & our.bm.legc['MPO',] > -9) | 
                                      (our.bm.legc['LYZ',] > -8.5 & our.bm.legc['HBB',] > -6) |
                                      (our.bm.legc['VPREB1',] > -12.5 & our.bm.legc['GATA1',] > -13.8) |
                                      (our.bm.legc['MPO',] > -8 & our.bm.legc['DNTT',] > -11)]
  # removing genes with very high ambience from feature list
  genes.for.proj.fil = setdiff(genes.for.proj, c('DNTT', 'MPO'))
  tmp.our.bm.cors = cor(bm.atlas.legc[genes.for.proj.fil,], our.bm.legc[genes.for.proj.fil,])
  our.bm.top5.matches = lapply(1:ncol(tmp.our.bm.cors), function(j) which(tmp.our.bm.cors[,j] >= sort(tmp.our.bm.cors[,j], decreasing=T)[5]))
  our.bm.proj.annotation = sapply(our.bm.top5.matches, function(x) most.common(bm.annotation[x]))
  names(our.bm.proj.annotation) = colnames(our.bm.legc)
  our.bm.proj.annotation[our.bm.legc['AVP',] + our.bm.legc['HLF',] > -28] = 'gold'
  our.bm.proj.annotation[our.bm.legc['SPIB',] > -13] = 'cyan'

  our.bm.projected.coords = do.call(rbind, lapply(1:length(our.bm.top5.matches), function(i) colMeans(umap.coords[colnames(bm.atlas.legc)[our.bm.top5.matches[[i]]],])))
  rownames(our.bm.projected.coords) = colnames(our.bm.legc)
  set.seed(43)
  umap.noise.x = rnorm(nrow(our.bm.projected.coords), 0, 0.3)
  umap.noise.y = rnorm(nrow(our.bm.projected.coords), 0, 0.3)
  names(umap.noise.x) = rownames(our.bm.projected.coords)
  names(umap.noise.y) = rownames(our.bm.projected.coords)
  #pal.mcs.to.project = colnames(pal.legc)[pal.legc['CD34',] > -15] #names(our.annotation)[our.annotation %in% c('plum', 'gold', 'steelblue', 'orange', 'chartreuse', 'blue', 'darkgreen')]
  our.bm.mcs.to.project = setdiff(colnames(our.bm.legc)[our.bm.legc['CD34',] > -14.5], our.bm.dbl.mcs)
  png(file.path(fig.dir, 'our_bm_and_bm_atlas_2d.png'), height=500, width=500)
  plot(bm.atlas.mdata$obs$umap_x, bm.atlas.mdata$obs$umap_y, bg='white', pch=21, cex=0.8, yaxt='n', xaxt='n')
  points(our.bm.projected.coords[our.bm.mcs.to.project,1] + umap.noise.x[our.bm.mcs.to.project], 
         our.bm.projected.coords[our.bm.mcs.to.project,2] + umap.noise.y[our.bm.mcs.to.project], bg=our.bm.proj.annotation[our.bm.mcs.to.project], pch=21, cex=1)
  dev.off()

  # gene scatters
  cur.fig.dir = file.path(fig.dir, 'bm_scatters')
  dir.create(cur.fig.dir, showWarnings=F)

  gene.pairs = list(c('MPO', 'AZU1'), c('DNTT', 'VPREB1'), c('GATA1', 'GATA2'), c('LMO4', 'HDC'), c('IRF8', 'SPIB'), c('MPO', 'DNTT'),
                    c('AVP', 'HOPX'), c('MME', 'JCHAIN'), c('RUNX3', 'HOXA9'), c('ZEB1', 'TCF7L2'))
  for (i in seq_along(gene.pairs)) {
    gene1 = gene.pairs[[i]][1]
    gene2 = gene.pairs[[i]][2]
    cur.xlim = range(our.legc[gene1, mcs.to.project], pal.legc[gene1, pal.mcs.to.project], our.bm.legc[gene1, our.bm.mcs.to.project])
    cur.ylim = range(our.legc[gene2, mcs.to.project], pal.legc[gene2, pal.mcs.to.project], our.bm.legc[gene2, our.bm.mcs.to.project])
    if (i < 7) {
      png(file.path(cur.fig.dir, sprintf('%s_%s.png', gene1, gene2)), height=500, width=1650)
      par(mfrow=c(1, 3))
    } else {
      png(file.path(cur.fig.dir, sprintf('%s_%s.png', gene1, gene2)), height=1650, width=500)
      par(mfrow=c(3, 1))
    }
    plot(our.legc[gene1, mcs.to.project], our.legc[gene2, mcs.to.project], xlab=gene1, ylab=gene2, xlim=cur.xlim, ylim=cur.ylim, cex=2.5, pch=21, bg=our.annotation[mcs.to.project])
    plot(our.bm.legc[gene1, our.bm.mcs.to.project], our.bm.legc[gene2, our.bm.mcs.to.project], xlab=gene1, ylab=gene2, xlim=cur.xlim, ylim=cur.ylim, cex=2.5, pch=21, bg=our.bm.proj.annotation[our.bm.mcs.to.project])
    plot(pal.legc[gene1, pal.mcs.to.project], pal.legc[gene2, pal.mcs.to.project], xlab=gene1, ylab=gene2, xlim=cur.xlim, ylim=cur.ylim, cex=2.5, pch=21, bg=pal.proj.annotation[pal.mcs.to.project])
    dev.off()
  }

  # HLF vs AVP
  #hlf.range = range(our.legc['HLF',], pbmc160.legc['HLF',], garvan.legc['HLF',], pal.legc['HLF',], bm.atlas.legc['HLF',])
  #avp.range = range(our.legc['AVP',], pbmc160.legc['AVP',], garvan.legc['AVP',], pal.legc['AVP',], bm.atlas.legc['AVP',])
  hlf.range = range(our.legc['HLF',], pal.legc['HLF',], bm.atlas.legc['HLF',], our.bm.legc['HLF',])
  avp.range = range(our.legc['AVP',], pal.legc['AVP',], bm.atlas.legc['AVP',], our.bm.legc['AVP',])
  #png(file.path(fig.dir, 'pbmc_bm_hsc.png'), height=500, width=2300)
  png(file.path(fig.dir, 'pbmc_bm_hsc.png'), height=500, width=2400)
  par(mfrow=c(1, 4))
  plot(our.legc['HLF',], our.legc['AVP',], xlim=hlf.range, ylim=avp.range, cex=3.5, pch=21, bg=our.annotation[colnames(our.legc)])
  plot(bm.atlas.legc['HLF',], bm.atlas.legc['AVP',], xlim=hlf.range, ylim=avp.range, cex=3.5, ylab='', pch=21, bg=bm.annotation[colnames(bm.atlas.legc)])
  plot(pal.legc['HLF',], pal.legc['AVP',], xlim=hlf.range, ylim=avp.range, cex=3.5, ylab='', pch=21, bg=pal.proj.annotation[colnames(pal.legc)])
  plot(our.bm.legc['HLF',], our.bm.legc['AVP',], xlim=hlf.range, ylim=avp.range, cex=3.5, ylab='', pch=21, bg=our.bm.proj.annotation[colnames(our.bm.legc)])
  #plot(1)
  #plot(garvan.legc['HLF',], garvan.legc['AVP',], pch=19, xlim=hlf.range, ylim=avp.range, cex=2.5)
  #plot(pbmc160.legc['HLF',], pbmc160.legc['AVP',], pch=19, xlim=hlf.range, ylim=avp.range, cex=2.5)
  dev.off()

}

bm.comparison.nktdp <- function(our.legc, fig.dir) {
  pbmc160.cdata = anndata::read_h5ad(file.path(MODEL.DIR, 'pbmc_160k_ref_cells.h5ad'))
  pbmc160.mdata = anndata::read_h5ad(file.path(MODEL.DIR, 'pbmc_160k_ref_metacells.h5ad'))
  pbmc160.legc = t(log2(1e-5 + pbmc160.mdata$X / rowSums(pbmc160.mdata$X)))

  garvan.cdata = anndata::read_h5ad(file.path(MODEL.DIR, 'pbmc_garvan_ref_cells.h5ad'))
  garvan.mdata = anndata::read_h5ad(file.path(MODEL.DIR, 'pbmc_garvan_ref_metacells.h5ad'))
  garvan.legc = t(log2(1e-5 + garvan.mdata$X / rowSums(garvan.mdata$X)))
  rownames(garvan.legc) = garvan.cdata$var$feature_name

  pal.cdata = anndata::read_h5ad(file.path(MODEL.DIR, 'bm_palantir_ref_cells.h5ad'))
  pal.mdata = anndata::read_h5ad(file.path(MODEL.DIR, 'bm_palantir_ref_metacells.h5ad'))
  pal.legc = t(log2(1e-5 + pal.mdata$X / rowSums(pal.mdata$X)))

  bm.atlas.mdata = anndata::read_h5ad(BM.MODEL.PATH)
  bm.atlas.legc = t(log2(1e-5 + bm.atlas.mdata$X / rowSums(bm.atlas.mdata$X)))

  our.bm.mdata = anndata::read_h5ad(file.path(MODEL.DIR, 'bm_our_ref_metacells.h5ad'))
  our.bm.mdata.x = our.bm.mdata$layers[['total_umis']]
  our.bm.legc = data.matrix(t(log2(1e-5 + our.bm.mdata.x / rowSums(our.bm.mdata.x))))
  colnames(our.bm.legc) = paste0('mc', 0:(ncol(our.bm.legc) - 1))

  acy3.range = range(our.legc['ACY3',], pbmc160.legc['ACY3',], garvan.legc['ACY3',], pal.legc['ACY3',], bm.atlas.legc['ACY3',], our.bm.legc['ACY3',])
  dntt.range = range(our.legc['DNTT',], pbmc160.legc['DNTT',], garvan.legc['DNTT',], pal.legc['DNTT',], bm.atlas.legc['DNTT',], our.bm.legc['DNTT',])
  syt2.vals = c(our.legc['SYT2',], pbmc160.legc['SYT2',], garvan.legc['SYT2',], pal.legc['SYT2',], bm.atlas.legc['SYT2',], our.bm.legc['SYT2',])
  syt2.range = range(syt2.vals)
  shades = colorRampPalette(c('darkblue', 'blue', 'white', 'red', 'darkred'))(100)
  png(file.path(fig.dir, 'pbmc_bm_acy3.png'), height=1100, width=1700)
  par(mfrow=c(2, 3))
  plot(our.bm.legc['ACY3',], our.bm.legc['DNTT',], pch=21, bg=shades[round(rescale(our.bm.legc['SYT2',], c(1, length(shades)), syt2.range))], xlim=acy3.range, ylim=dntt.range, cex=3.5)
  plot(pal.legc['ACY3',], pal.legc['DNTT',], pch=21, bg=shades[round(rescale(pal.legc['SYT2',], c(1, length(shades)), syt2.range))], xlim=acy3.range, ylim=dntt.range, cex=3.5)
  plot(bm.atlas.legc['ACY3',], bm.atlas.legc['DNTT',], pch=21, bg=shades[round(rescale(bm.atlas.legc['SYT2',], c(1, length(shades)), syt2.range))], xlim=acy3.range, ylim=dntt.range, cex=3.5)
  plot(our.legc['ACY3',], our.legc['DNTT',], pch=21, bg=shades[round(rescale(our.legc['SYT2',], c(1, length(shades)), syt2.range))], xlim=acy3.range, ylim=dntt.range, cex=3.5)
  plot(garvan.legc['ACY3',], garvan.legc['DNTT',], pch=21, bg=shades[round(rescale(garvan.legc['SYT2',], c(1, length(shades)), syt2.range))], xlim=acy3.range, ylim=dntt.range, cex=3.5)
  plot(pbmc160.legc['ACY3',], pbmc160.legc['DNTT',], pch=21, bg=shades[round(rescale(pbmc160.legc['SYT2',], c(1, length(shades)), syt2.range))], xlim=acy3.range, ylim=dntt.range, cex=3.5)
  dev.off()
  # plot legend
  #image.plot(matrix(syt2.vals, nrow=2), col=shades)
}

bemp.nktdp.scatters <- function(our.legc, our.annotation, fig.dir=BASE.FIG.DIR) {
  
  cur.mcs = names(our.annotation[our.annotation%in% c('purple', 'plum')])
  genes.for.bemp.scatters = c('GATA2', 'CREB3L2', 'LMO2', 'POU2F2', 'TLE4', 'KLF1', 
                              'CNRIP1', 'HPGDS', 'TET2', 'TNFSF10', 'CD34', 'HBD', 'CD74', 'BLVRB')
  cur.fig.dir = file.path(fig.dir, 'bemp_figs')
  dir.create(cur.fig.dir, showWarnings=F)
  for(g in genes.for.bemp.scatters) {
    png(file.path(cur.fig.dir, paste0(g, '.png')), w=500 ,h=300)
    par(mar=c(0,4,0,1))
    plot(our.legc['LMO4', cur.mcs], our.legc[g, cur.mcs], pch=19, col=our.annotation[cur.mcs], xlab='LMO4', ylab=g, cex=2)
    dev.off()
  }

  f = our.annotation == 'lightblue' & our.legc['IRF8',] + our.legc['TCF7',] > -28.7
  grad = our.legc["TCF7",f]-our.legc["IRF8",f]
  c_tcf7irf8 =  apply(our.legc[,f],1,cor,grad,m="spearman")

  cur.fig.dir = file.path(fig.dir, 'acy3_grad_figs')
  dir.create(cur.fig.dir, showWarnings=F)
  genes.for.nktdp.scatters = c('CD7', 'MAF', 'SPI1', 'TRBC2', 'IL7R', 'CD74',
                               'HLA-DMA', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DQB1', 'HLA-DRA', 'HLA-DRB1')

  for(g in genes.for.nktdp.scatters) {
    png(file.path(cur.fig.dir, paste0(g, '.png')), w=500,h=300)
    par(mar=c(0,4,0,1))
    plot(grad, our.legc[g,f], cex=2, pch=19, col="lightblue", ylab=g)
    dev.off()
  }


}

plot.hsc.modes.and.dynamics <- function(our.legc, our.annotation, fig.dir=BASE.FIG.DIR) {
  col_clp5 = colorRampPalette(c("brown", "blue"))(5)
  col_mep4 = colorRampPalette(c("gold", "plum"))(4)

  clp.mcs.for.modes = colnames(our.legc)[our.annotation %in% c('blue', 'darkblue', '#6E1C71') & our.legc['AVP',] > -15 & our.legc['AVP',] < -13]
  mep.mcs.for.modes = colnames(our.legc)[our.annotation %in% c('gold', '#EEBB6E', '#plum', '##7F9D00') & our.legc['AVP',] > -15 & our.legc['AVP',] < -13]
  e_clp_for_modes = rowMeans(our.legc[,clp.mcs.for.modes]) 
  e_mep_for_modes = rowMeans(our.legc[,mep.mcs.for.modes]) 

  e_hsc = rowMeans(our.legc[,our.annotation[colnames(our.legc)] == 'brown'])
  e_mep_e = rowMeans(our.legc[,our.annotation[colnames(our.legc)] == '#EEBB6E'])
  e_clp_multi = rowMeans(our.legc[,our.annotation[colnames(our.legc)] == 'blue'])

  #dlt_hsc_clp = e_hsc - e_clp_for_modes
  #dlt_hsc_mep = e_hsc - e_mep_for_modes
  dlt_hsc_clp = e_hsc - e_clp_multi
  dlt_hsc_mep = e_hsc - e_mep_e

  diff_thresh = log2(1.75)
  hsc_mode = ifelse(dlt_hsc_clp > diff_thresh & dlt_hsc_mep > diff_thresh, "I",
                                        ifelse(dlt_hsc_clp < -diff_thresh & dlt_hsc_mep < -diff_thresh, "II",
                                        ifelse(dlt_hsc_clp >  diff_thresh & dlt_hsc_mep < -diff_thresh, "IV",
                                        ifelse(dlt_hsc_clp < -diff_thresh & dlt_hsc_mep >  diff_thresh, "III",
                                        ifelse(dlt_hsc_clp < -diff_thresh, "V",
                                        ifelse(dlt_hsc_mep < -diff_thresh, "VI",
                                        ifelse(dlt_hsc_clp >  diff_thresh, "VII",
                                        ifelse(dlt_hsc_mep >  diff_thresh, "VIII","0"))))))))

  print(table(hsc_mode))

  tmp_lim = range(dlt_hsc_mep, dlt_hsc_clp)
  png(file.path(fig.dir, 'hsc_genemap.png'), w=400, h=400)
  plot(dlt_hsc_mep, dlt_hsc_clp, pch=19, xlab=NA, ylab=NA, xlim=tmp_lim, ylim=tmp_lim,
                          col=ifelse(dlt_hsc_clp >  diff_thresh & dlt_hsc_mep >  diff_thresh, "brown",
                              ifelse(dlt_hsc_clp < -diff_thresh & dlt_hsc_mep < -diff_thresh, "brown",
                              ifelse(dlt_hsc_clp >  diff_thresh & dlt_hsc_mep < -diff_thresh, col_mep4[4],
                              ifelse(dlt_hsc_clp < -diff_thresh & dlt_hsc_mep >  diff_thresh, col_clp5[5],
                              ifelse(dlt_hsc_clp < -diff_thresh, col_clp5[3],
                              ifelse(dlt_hsc_mep < -diff_thresh, col_mep4[3], "gray")))))))
  abline(h= diff_thresh, lty=2)
  abline(h=-diff_thresh, lty=2)
  abline(v= diff_thresh, lty=2)
  abline(v=-diff_thresh, lty=2)
  dev.off()

  hsc.diff.tbl = data.frame(gene=names(hsc_mode), hsc_exp=e_hsc, clpm_exp=e_clp_multi, mebempe_exp=e_mep_e, hsc_clpm_delta=dlt_hsc_clp, hsc_mebempe_delta=dlt_hsc_mep, class=hsc_mode)
  write.csv(hsc.diff.tbl, file=file.path(SUPP.TABLE.DIR, 's4_hsc_expression_diffs.csv'), quote=F, row.names=F)

  tfs = read.table(TFS.FILE.PATH, stringsAsFactors=F)
  tfs = toupper(tfs$x)
  tfs = intersect(tfs, rownames(our.legc))
  print(sort(intersect(names(which(dlt_hsc_mep > log2(1.25) & dlt_hsc_clp > log2(1.25))), tfs)))

  # and now the trajectories
  f_filt_pls = ((our.annotation %in% c('brown', 'gold', '#EEBB6E', 'plum', 'red', '#7F9D00')) & (our.legc['HOXA9',] < -13.1) | 
                (our.annotation %in% c('#6E1C71', 'blue', 'darkblue')) & (our.legc['HOXA9',] >= -13.1))
  
  traj_ord = order(ifelse(our.legc["HOXA9",f_filt_pls] > -13.1,10 - our.legc["AVP",f_filt_pls],our.legc["AVP",f_filt_pls]))
  n_mep = sum(our.legc["HOXA9",f_filt_pls]< -13.1)
  n_clp = sum(our.legc["HOXA9",f_filt_pls]> -13.1)
  
  plt_traj_gene = function(g, n_window=30)
  {
          y = our.legc[g,f_filt_pls][traj_ord]
          sy = rollmean(y, n_window, fill='extend')
          x_traj = c(seq(0,2,l=n_mep),seq(2,3, l=n_clp))
          yl = c(pmin(pmax(e_hsc[g]-2,log2(0.5e-5)),min(y)),
                          pmax(e_hsc[g]+2,max(y)))
          plot(x_traj, y , pch=19, col=our.annotation[f_filt_pls][traj_ord], cex=0.8, xlab=NA, ylab=g, ylim=yl)
          lines(x_traj, sy, pch=19, col="black", cex=0.2, lwd=2)
          abline(h=e_hsc[g])
          abline(h=e_hsc[g]-1, lty=2)
          abline(h=e_hsc[g]+1, lty=2)
  }
  
  high_hsc_tfs = names(which(e_hsc[tfs]> -13))
  diff_grad_tfs = names(which(abs(e_mep_e-e_clp_multi)[tfs]>1))
  diff_hsc_tfs = union(names(which(abs(e_hsc-e_clp_multi)[tfs]>0.9)),
                       names(which(abs(e_hsc-e_mep_e)[tfs]>0.8)))
  
  cand_tfs = union(high_hsc_tfs, union(diff_grad_tfs, diff_hsc_tfs))
  
  tfs.fig.dir = file.path(fig.dir, 'grad_tfs')
  dir.create(tfs.fig.dir, showWarnings=F)
  for(tf in cand_tfs) {
    png(sprintf(file.path(tfs.fig.dir, '%s.png'), tf), w=800, h=300)
    par(mar=c(0, 4, 0, 1))
    plt_traj_gene(tf)
    dev.off()
  }
}

get.age.and.sex <- function(our.cdata, indiv.to.used.samp, cell.to.samp) {
  stopifnot(length(cell.to.samp) == nrow(our.cdata$obs))
  stopifnot(all(names(cell.to.samp) == rownames(our.cdata$obs)))

  samp.src.tbl = table(cell.to.samp, as.character(our.cdata$obs$src))
  stopifnot(all(rowSums(samp.src.tbl > 0) == 1))
  samp.to.src = colnames(samp.src.tbl)[apply(samp.src.tbl, 1, which.max)]
  names(samp.to.src) = rownames(samp.src.tbl)

  age.and.sex = read.csv(ORIG.SEX.AGE.PATH, stringsAsFactors=F)
  age.and.sex = age.and.sex[!is.na(age.and.sex$num),]
  rownames(age.and.sex) = age.and.sex$num
  sexes = age.and.sex[,'gender']
  ages = age.and.sex[,'age']
  names(sexes) = paste0('N', rownames(age.and.sex))
  names(ages) = paste0('N', rownames(age.and.sex))

  new.age.and.sex = read.csv(NEW.SEX.AGE.PATH, header=T)
  rownames(new.age.and.sex) = new.age.and.sex$num
  new.sexes = new.age.and.sex[,'gender']
  new.ages = new.age.and.sex[,'age']
  names(new.sexes) = rownames(new.age.and.sex)
  names(new.ages) = rownames(new.age.and.sex)

  common.sex.indivs = intersect(names(sexes), names(new.sexes))
  stopifnot(all(sexes[common.sex.indivs] == new.sexes[common.sex.indivs]))
  new.sexes.fil = new.sexes[setdiff(names(new.sexes), names(sexes))]
  returned.sexes = c(sexes, new.sexes.fil)
  stopifnot(all(names(indiv.to.used.samp) %in% names(returned.sexes)))


  exception.indivs = c('N171', 'N180', 'N91', 'N133', 'N130', 'N135', 'N152', 'N153', 'N156', 'N155', 'N159', 'N167', 'N169', 'N168', 'N176', 'N166')
  indivs.with.old.samp = c(names(indiv.to.used.samp)[samp.to.src[indiv.to.used.samp] != 'new'], exception.indivs)
  indivs.with.new.samp = setdiff(names(indiv.to.used.samp)[samp.to.src[indiv.to.used.samp] == 'new'], exception.indivs)
  stopifnot(length(intersect(indivs.with.old.samp, indivs.with.new.samp)) == 0)

  stopifnot(all(indivs.with.old.samp %in% names(ages)))
  stopifnot(all(indivs.with.new.samp %in% names(new.ages)))
  returned.ages = c(ages[indivs.with.old.samp], new.ages[indivs.with.new.samp])

  return(list(ages=returned.ages[names(indiv.to.used.samp)], sexes=returned.sexes[names(indiv.to.used.samp)]))
}

mk.model.analysis <- function(fig.dir=file.path(BASE.FIG.DIR, 'fig1')) {
  mdata.output.all = anndata::read_h5ad(file.path(MODEL.DIR, '148_indiv_ref_metacells_unfiltered.h5ad'))
  tmp.umi.counts = mdata.output.all$layers[['total_umis']]
  doublet.legc = log2(1e-5 + t(tmp.umi.counts / rowSums(tmp.umi.counts)))

  # doublet heatmap
  genes.for.heatmap = c('PF4', 'PPBP', 'GATA1', 'LYZ', 'DNTT', 'CD3D', 'CD79A')
  tmp.legc = doublet.legc[genes.for.heatmap,]
  sus.mk.mcs = colnames(tmp.legc)[tmp.legc['PF4',] > -11.5 | tmp.legc['PPBP',] > -10.2]
  tmp.legc.norm = tmp.legc - rowMeans(tmp.legc)
  shades = colorRampPalette(c("darkblue", "blue", "white",'red', "darkred"))(101)
  pheatmap(t(pmin(pmax(tmp.legc.norm[,sus.mk.mcs], -2), 2)), col=shades, 
           breaks=seq(-2, 2, length.out=length(shades) + 1), cluster_cols=F, filename=file.path(fig.dir, 'mk_doublets.png'), height=5, width=5, show_rownames=F)
}

most.common <- function(x) {
  tmp.tbl = table(x)
  return(names(tmp.tbl)[which.max(tmp.tbl)])
}

plot.full.umap <- function(our.legc, our.annotation, fig.dir) {
  remove_types = c("white", "grey")
  mc_col = our.annotation
  f_filt = !mc_col %in% remove_types
  
  c_vim = apply(our.legc[,f_filt],1,cor,our.legc["VIM", f_filt])
  
  g_var_stem = apply(our.legc[,f_filt], 1, function(x) max(x)-min(x))
  
  g_feat_stemp = names(which(g_var_stem>2 & apply(our.legc[,f_filt],1,max)>-14 & abs(c_vim)<0.33))
  
  lfp = our.legc - rowMeans(our.legc)
  used.ctypes = setdiff(unique(our.annotation), remove_types)
  g_feat_stemp = unique(unlist(lapply(used.ctypes, function(cur.ctype) {
    cur.feats = names(tail(sort(rowMeans(lfp[,mc_col==cur.ctype])),30))
    return(cur.feats)
  })))
  
  g_feat_stemp = setdiff(g_feat_stemp, grep("^RPL", g_feat_stemp,v=T))
  g_feat_stemp = setdiff(g_feat_stemp, grep("^RPS", g_feat_stemp,v=T))
  g_feat_stemp = setdiff(g_feat_stemp, grep("^MRPL", g_feat_stemp,v=T))
  g_feat_stemp = setdiff(g_feat_stemp, c("PRSS2"))
  g_feat_stemp = g_feat_stemp[c_vim[g_feat_stemp] < 0.4]
  
  conf = umap::umap.defaults
  conf$min_dist=0.96
  conf$random_state = 45
  conf$n_neighbors=6
  conf$n_epoch=2000
  
  uf = t(our.legc[g_feat_stemp, f_filt])
  uf = pmax(uf, -14)
  um = umap::umap(uf, config=conf)
  
  plot(um$layout[,1], um$layout[,2], col=mc_col[f_filt], pch=19,cex=1.5)
  
  png(file.path(fig.dir, "full_umap.png"), w=800, h=800)
  plot(um$layout[,1], um$layout[,2], col=mc_col[f_filt], pch=19,cex=1.5, yaxt='n', xaxt='n')
  dev.off()

  return(um$layout)

}

plot.filtered.umap <- function(our.legc, our.annotation, fig.dir) {
  remove_types = c("orange", "cyan","chartreuse", "steelblue", "darkgreen", "black", "white", "grey")
  mc_col = our.annotation
  f_filt = !mc_col %in% remove_types
  
  c_vim = apply(our.legc[,f_filt],1,cor,our.legc["VIM", f_filt])
  
  g_var_stem = apply(our.legc[,f_filt], 1, function(x) max(x)-min(x))
  
  g_feat_stemp = names(which(g_var_stem>2 & apply(our.legc[,f_filt],1,max)>-14 & abs(c_vim)<0.33))
  
  lfp = our.legc - rowMeans(our.legc)
  
  hsc_feat = names(tail(sort(rowMeans(lfp[,mc_col=="brown"])),30))
  b_feat = names(tail(sort(rowMeans(lfp[,mc_col=="darkblue"])),30))
  #acy_feat = names(tail(sort(rowMeans(lfp[,mc_col=="lightblue"])),30))
  acy_feat = names(tail(sort(rowMeans(lfp[,mc_col=="lightblue"])),50))
  #dntt_feat = names(tail(sort(rowMeans(lfp[,mc_col=="blue"])),30))
  dntt_feat = names(tail(sort(rowMeans(lfp[,mc_col=="blue"])),50))
  mep_feat = names(tail(sort(rowMeans(lfp[,mc_col=="plum"])),30))
  mpo_feat = names(tail(sort(rowMeans(lfp[,mc_col=="darkgreen"])),30))
  baso_feat = names(tail(sort(rowMeans(lfp[,mc_col=="purple"])),30))
  
  g_feat_stemp = unique(c(hsc_feat, b_feat, acy_feat, mep_feat, mpo_feat, baso_feat, dntt_feat))
  ery_feat = names(tail(sort(rowMeans(lfp[,mc_col=="red"])),30))
  g_feat_stemp = c(g_feat_stemp, ery_feat)
  
  g_feat_stemp = setdiff(g_feat_stemp, grep("^RPL", g_feat_stemp,v=T))
  g_feat_stemp = setdiff(g_feat_stemp, grep("^RPS", g_feat_stemp,v=T))
  g_feat_stemp = setdiff(g_feat_stemp, grep("^MRPL", g_feat_stemp,v=T))
  g_feat_stemp = setdiff(g_feat_stemp, c("PRSS2"))
  g_feat_stemp = g_feat_stemp[c_vim[g_feat_stemp] < 0.4]
  
  conf = umap::umap.defaults
  conf$min_dist=0.96
  conf$random_state = 45
  conf$n_neighbors=6
  conf$n_epoch=2000
  
  f_filt[our.legc["VPREB3",]> -13]= F
  f_filt[our.legc["HBB",]> -10]= F
  uf = t(our.legc[g_feat_stemp, f_filt])
  uf = pmax(uf, -14)
  um = umap::umap(uf, config=conf)
  
  if (is.null(fig.dir)) {
    plot(um$layout[,1], um$layout[,2], col=mc_col[f_filt], pch=19,cex=1.5)
  } else {
    png(file.path(fig.dir, "filtered_umap.png"), w=600, h=800)
    plot(um$layout[,1], um$layout[,2], col=mc_col[f_filt], pch=19,cex=1.5)
    dev.off()
  }

  return(um$layout)

}

rescale <- function(x, dest.range, src.range=NULL) {
  from = dest.range[1]
  to = dest.range[2]
  if (is.null(src.range)) {
    maxx = max(x)
    minx = min(x)
  } else {
    maxx = src.range[2]
    minx = src.range[1]
  }
  out = (to - from) * (x - minx)
  out = out / (maxx - minx)
  out + from
}

annotate.model <- function(our.legc) {

  annotation.thresholds = list(
                            list(gene1='ADIRF', gene2='CAV1', thresh1=-14, thresh2=-20, logical.operator='and', cur.col='black'),
                            list(gene1='IRF8', gene2='SPIB', thresh1=-11.3, thresh2=-13, logical.operator='and', cur.col='cyan'),
                            list(gene1='NKG7', gene2='CD3D', thresh1=-9, thresh2=-13, logical.operator='or', cur.col='orange'),
                            list(gene1='TCF7', gene2='CD3D', thresh1=-12, thresh2=-20, logical.operator='and', cur.col='orange'),
                            list(gene1='AVP', gene2='HLF', thresh1=-25, thresh2=NA, logical.operator='sum', cur.col='brown'),
                            list(gene1='AVP', gene2='HLF', thresh1=-10.6, thresh2=-20, logical.operator='and', cur.col='brown'),
                            list(gene1='S100A8', gene2='S100A9', thresh1=-13, thresh2=-20, logical.operator='and', cur.col='chartreuse'),
                            list(gene1='LMO4', gene2='HDC', thresh1=-22, thresh2=NA, logical.operator='sum', cur.col='purple'),
                            list(gene1='HBB', gene2='HBD', thresh1=-11.5, thresh2=-10.3, logical.operator='or', cur.col='red'),
                            list(gene1='HBB', gene2='HBD', thresh1=-26.5, thresh2=NA, logical.operator='sum', cur.col='red'),
                            list(gene1='MPO', gene2='AZU1', thresh1=-10, thresh2=-20, logical.operator='and', cur.col='darkgreen')
			  )

  mc.cols = rep(NA, ncol(our.legc))
  names(mc.cols) = colnames(our.legc)
  for (i in seq_along(annotation.thresholds)) {
    gene1 = annotation.thresholds[[i]]$gene1
    gene2 = annotation.thresholds[[i]]$gene2
    thresh1 = annotation.thresholds[[i]]$thresh1
    thresh2 = annotation.thresholds[[i]]$thresh2
    cur.logical = annotation.thresholds[[i]]$logical.operator
    cur.col = annotation.thresholds[[i]]$cur.col
    
    #plot(our.legc[gene1,], our.legc[gene2,], pch=21, bg=our.annotation, xlab=gene1, ylab=gene2)
    #abline(v=thresh1, col=2)
    #abline(h=thresh2, col=2)

    if (cur.logical == 'and') {
      cur.mcs = colnames(our.legc)[our.legc[gene1,] > thresh1 & our.legc[gene2,] > thresh2]
    } else if (cur.logical == 'or') {
      cur.mcs = colnames(our.legc)[our.legc[gene1,] > thresh1 | our.legc[gene2,] > thresh2]
    } else if (cur.logical == 'sum') {
      cur.mcs = colnames(our.legc)[our.legc[gene1,] + our.legc[gene2,] > thresh1]
    }
    cur.na.mcs = names(mc.cols)[is.na(mc.cols)]
    #stopifnot(all(cur.mcs %in% cur.na.mcs))

    #plot(our.legc[gene1,], our.legc[gene2,], pch=21, bg=(colnames(our.legc) %in% cur.mcs) + 1, xlab=gene1, ylab=gene2)

    mc.cols[intersect(cur.na.mcs, cur.mcs)] = cur.col
  }

  mc.cols[our.legc['MPO',] > -13.5 & !(mc.cols %in% c('darkgreen'))] = '#7F9D00'
  mc.cols[our.legc['CD79A',] > -12 & our.legc['VPREB1',] < -12] = 'steelblue'
  mc.cols[our.legc['RUNX3',] + our.legc['ACY3',] > -25 & is.na(mc.cols)] = 'lightblue'
  mc.cols[our.legc['HOXA9',] + our.legc['QPRT',] > -27.3 & !(mc.cols %in% c('brown', 'lightblue'))] = '#6E1C71'
  mc.cols[our.legc['HOXA9',] + our.legc['QPRT',] > -25 & !(mc.cols %in% c('brown', 'lightblue', 'darkblue'))] = 'blue'
  #mc.cols[our.legc['DNTT',] + our.legc['VPREB1',] > -25 & !(mc.cols %in% c('#6E1C71'))] = 'darkblue'
  mc.cols[our.legc['DNTT',] + our.legc['VPREB1',] > -24] = 'darkblue'
  mc.cols[our.legc['GATA1',] + our.legc['KLF1',] > -26.3 & is.na(mc.cols)] = 'plum'
  mc.cols[our.legc['GATA1',] + our.legc['KLF1',] > -30.7 & is.na(mc.cols)] = '#EEBB6E'

  gata2.cors = apply(our.legc, 1, cor, our.legc['GATA2',])
  hoxa9.cors = apply(our.legc, 1, cor, our.legc['HOXA9',])
  gata2.genes = names(tail(sort(gata2.cors), n=50))
  hoxa9.genes = names(tail(sort(hoxa9.cors), n=50))
  gata2.scores = colMeans(our.legc[gata2.genes,])
  hoxa9.scores = colMeans(our.legc[hoxa9.genes,])

  stopifnot(all(!is.na(mc.cols[gata2.scores < -14.2 & hoxa9.scores > -14])))
  mc.cols[gata2.scores > -14.2 & hoxa9.scores > -14 & is.na(mc.cols)] = '#6E1C71'
  mc.cols[gata2.scores > -14.5 & is.na(mc.cols)] = 'gold'
  mc.cols[our.legc['CD79A',] > -12 & is.na(mc.cols)] = 'steelblue'
  mc.cols[our.legc['LYZ',] > -10 & is.na(mc.cols)] = 'chartreuse'

  mc.cols[is.na(mc.cols)] = 'grey'
  #plot(gata2.scores, hoxa9.scores, pch=21, bg=mc.cols)
  mc.cols[our.legc['CD79A',] > -13.5 & mc.cols %in% c('chartreuse', '#EEBB6E', '#6E1C71')] = 'grey'
  mc.cols[our.legc['ACY3',] > -13 & mc.cols %in% c('chartreuse', 'gold')] = 'grey'
  mc.cols[our.legc['GATA1',] > -14.5 & mc.cols %in% c('chartreuse')] = 'grey'

  mc.cols[our.legc['HOXA9',] > -13 & our.legc['GATA1',] < -15 & mc.cols %in% c('#EEBB6E', 'gold')] = '#6E1C71'

  # and now doublet detection - grey is doublet from now on
  doublet.thresholds = list(
                            list(gene1='LYZ', gene2='HBB', thresh1=-13.5, thresh2=-12.5),
                            list(gene1='LYZ', gene2='HBD', thresh1=-10, thresh2=-14),
                            list(gene1='LYZ', gene2='DNTT', thresh1=-12, thresh2=-13),
                            list(gene1='CD3D', gene2='LYZ', thresh1=-14.5, thresh2=-12),
                            list(gene1='GATA1', gene2='LYZ', thresh1=-14.7, thresh2=-12),
                            list(gene1='GATA1', gene2='ACY3', thresh1=-14, thresh2=-12),
                            list(gene1='GATA2', gene2='CD79A', thresh1=-14, thresh2=-13),
                            list(gene1='AVP', gene2='LYZ', thresh1=-13, thresh2=-12),
                            list(gene1='AVP', gene2='CD79A', thresh1=-13, thresh2=-13.5),
                            list(gene1='AVP', gene2='ACY3', thresh1=-13, thresh2=-12.5),
                            list(gene1='CD3D', gene2='LYZ', thresh1=-13, thresh2=-12),
                            list(gene1='CD3D', gene2='CD79A', thresh1=-14, thresh2=-13),
                            list(gene1='CD3D', gene2='ADIRF', thresh1=-14, thresh2=-12),
                            list(gene1='MPO', gene2='CD79A', thresh1=-13, thresh2=-13),
                            list(gene1='S100A8', gene2='CD79A', thresh1=-11, thresh2=-13),
                            list(gene1='MPO', gene2='RUNX3', thresh1=-13.8, thresh2=-13.3),
                            list(gene1='MPO', gene2='RUNX3', thresh1=-13, thresh2=-14),
                            list(gene1='HBD', gene2='RUNX3', thresh1=-13.5, thresh2=-13),
                            list(gene1='MPO', gene2='ACY3', thresh1=-12, thresh2=-12),
                            list(gene1='CD79A', gene2='LYZ', thresh1=-12, thresh2=-12),
                            list(gene1='GATA1', gene2='LYZ', thresh1=-14.5, thresh2=-12.5)
			   )
  doublet.mcs = c()
  already.doublets = names(mc.cols)[mc.cols == 'grey']
  #par(mfrow=c(1, 3))
  for (i in seq_along(doublet.thresholds)) {
  #for (i in poss.doublets) {
  #for (i in 13:18) {
    #par(mfrow=c(1, 2))
    #par(mfrow=c(1, 3))
    cur.thresh = doublet.thresholds[[i]]
    #cur.thresh = doublet.thresholds.full[[i]]
    gene1 = cur.thresh$gene1
    gene2 = cur.thresh$gene2
    thresh1 = cur.thresh$thresh1
    thresh2 = cur.thresh$thresh2
    cur.dbl.mcs = colnames(our.legc)[our.legc[gene1, ] > thresh1 & our.legc[gene2, ] > thresh2]
    doublet.mcs = c(doublet.mcs, cur.dbl.mcs)
    next

    cur.xlim = range(our.legc.prev[gene1,], our.legc[gene1,])
    cur.ylim = range(our.legc.prev[gene2,], our.legc[gene2,])
    plot(our.legc.prev[gene1,], our.legc.prev[gene2,], bg=our.annotation.prev, pch=21, xlab=gene1, ylab=gene2, main='prev model', xlim=cur.xlim, ylim=cur.ylim)
    abline(v=thresh1)
    abline(h=thresh2)
    plot(our.legc[gene1,], our.legc[gene2,], bg=(colnames(our.legc) %in% c(doublet.mcs, already.doublets)) + 1, pch=21, xlab=gene1, ylab=gene2, main='current model', xlim=cur.xlim, ylim=cur.ylim)
    abline(v=thresh1)
    abline(h=thresh2)
    plot(our.legc[gene1,], our.legc[gene2,], bg=mc.cols, pch=21, xlab=gene1, ylab=gene2, main='current model', xlim=cur.xlim, ylim=cur.ylim)
    abline(v=thresh1)
    abline(h=thresh2)
  }
  mc.cols[doublet.mcs] = 'grey'

  return(mc.cols)
}

