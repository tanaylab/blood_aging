library(glmnet)

fig2 <- function() {
  # exactly like the beginning of fig1...
  dir.create(BASE.FIG.DIR, showWarnings=F)
  fig.dir = file.path(BASE.FIG.DIR, 'fig2')
  dir.create(fig.dir, showWarnings=F)

  our.cdata = anndata::read_h5ad(file.path(MODEL.DIR, '148_indiv_ref_cells.h5ad'))
  our.mdata = anndata::read_h5ad(file.path(MODEL.DIR, '148_indiv_ref_metacells.h5ad'))
  tmp.umi.counts = our.mdata$layers[['total_umis']]
  our.legc = log2(1e-5 + t(tmp.umi.counts / rowSums(tmp.umi.counts)))
  colnames(our.legc) = paste0('mc', 0:(ncol(our.legc) - 1))
  our.assign = paste0('mc', our.cdata$obs$metacell)
  names(our.assign) = rownames(our.cdata$obs)

  our.mdata.x = our.mdata$layers[['total_umis']]
  rownames(our.mdata.x) = paste0('mc', 0:(nrow(our.mdata.x) - 1))

  our.annotation = annotate.model(our.legc)

  # distribution across metacells
  cell.to.meas = sprintf('%s-%s', our.cdata$obs$exp_name, our.cdata$obs$indiv_id)
  names(cell.to.meas) = rownames(our.cdata$obs)
  cell.to.samp = gsub('_\\d_ultima-', '_ultima-', gsub('_\\d-', '-', cell.to.meas))
  tmp.tbl = table(cell.to.samp, our.annotation[as.character(our.assign[names(cell.to.samp)])])
  all.samp = unique(cell.to.samp)
  indiv.to.used.samp = get.one.samp.per.indiv(all.samp, our.cdata, cell.to.samp)
  samp.to.indiv = sapply(strsplit(all.samp, '-'), function(x) x[2])
  names(samp.to.indiv) = all.samp

  age.and.sex = get.age.and.sex(our.cdata, indiv.to.used.samp, cell.to.samp)
  ages = age.and.sex$ages
  sexes = age.and.sex$sexes

  all.indivs = as.character(unique(our.cdata$obs$indiv_id))
  is.male = sexes == 'male'
  all.males = all.indivs[sexes[all.indivs] == 'male']
  all.females = all.indivs[sexes[all.indivs] == 'female']
  indivs.cols = ifelse(is.male[all.indivs], 'lightgreen', '#FFD580')

  # edf  4A, 4B
  de.between.indivs.analysis(our.cdata, our.mdata, our.annotation, our.assign, indiv.to.used.samp, cell.to.samp, fig.dir=fig.dir)

  # individuals for which cellranger cut some cells with low number of UMIs, potentially lowering their fraction of CLPs
  # (which have a lower number of UMIs)
  indiv.to.used.samp.fil = indiv.to.used.samp[!(names(indiv.to.used.samp) %in% c('N288', 'N290', 'N365'))]

  # fig 2B
  plot.cell.type.freqs(our.assign, our.annotation, our.legc, cell.to.samp, indiv.to.used.samp.fil, fig.dir=fig.dir)

  # main trajectory analysis
  # fig 2D
  plot.main.traj(our.legc, our.annotation, cell.to.samp, indiv.to.used.samp.fil, ages, indivs.cols)

  # biological and technical replicates
  # illumina - ultima technical similarity
  tech.rep.pools = get.tech.rep.pools(our.cdata$X, our.cdata$obs, rownames(our.cdata$obs))
  rep.illum.pool = tech.rep.pools$rep.illum.pool
  rep.ultima.pool = tech.rep.pools$rep.ultima.pool
  rep.illum.legc = log2(1e-5 + t(t(rep.illum.pool) / colSums(rep.illum.pool)))
  rep.ultima.legc = log2(1e-5 + t(t(rep.ultima.pool) / colSums(rep.ultima.pool)))

  # edf 7A
  cur.lim = range(c(rep.illum.legc, rep.ultima.legc))
  png(file.path(fig.dir, 'ultima_qc.png'), width=3800, height=950)
  par(mfrow=c(2, 6))
  for (i in 1:12) {
    #plot(rep.illum.legc[,i], rep.ultima.legc[,i], pch=19, cex=0.8, xlim=cur.lim, ylim=cur.lim)
    plot(rep.illum.legc[,i], rep.ultima.legc[,i], pch=19, cex=1.5, xlim=cur.lim, ylim=cur.lim)
  }
  dev.off()

  # technical and biological replicates
  tech.reps = get.tech.replicates(cell.to.meas, names(cell.to.meas))
  bio.reps = get.biological.replicates(all.samp, indiv.to.used.samp.fil)
  ctype.to.col = c("CLP-E"="#6E1C71", "GMP-E"="#7F9D00", "MEBEMP-E"="#EEBB6E", "Endothel"="black", "CLP-M"="blue", "HSC"="brown",
                   "Monocytes"="chartreuse", "DC"="cyan", "CLP-L"="darkblue", "GMP-L"="darkgreen", "MPP"="gold", "NKTDP"="lightblue", "NKT"="orange",
                   "MEBEMP-L"="plum", "BEMP"="purple", "EP"="red", "B"="steelblue", "Unknown"="white")
  col.to.ctype = names(ctype.to.col)
  names(col.to.ctype) = ctype.to.col

  cur.used.ctypes = c('purple', 'red', 'plum', '#EEBB6E', "#7F9D00", 'gold', 'brown', "#6E1C71", 'blue', 'darkblue', 'lightblue')
  #samp.ctype.tbl = get.ctype.distrib(our.assign, our.annotation, cell.to.samp)[,cur.used.ctypes]
  #meas.ctype.tbl = get.ctype.distrib(our.assign, our.annotation, cell.to.meas)[,cur.used.ctypes]
  cd34.pos.mcs = colnames(our.legc)[our.legc['CD34',] > -14.5]
  samp.ctype.tbl = get.ctype.distrib(our.assign, our.annotation, cell.to.samp, cd34.pos.mcs)[,cur.used.ctypes]
  meas.ctype.tbl = get.ctype.distrib(our.assign, our.annotation, cell.to.meas, cd34.pos.mcs)[,cur.used.ctypes]
  should.use.tech.rep = rowSums(meas.ctype.tbl[unlist(tech.reps[,1]),]) >= 500 & rowSums(meas.ctype.tbl[unlist(tech.reps[,2]),]) >= 500
  should.use.bio.rep = rowSums(samp.ctype.tbl[unlist(bio.reps[,1]),]) >= 500 & rowSums(samp.ctype.tbl[unlist(bio.reps[,2]),]) >= 500
  samp.ctype.tbl.norm = samp.ctype.tbl / rowSums(samp.ctype.tbl)
  meas.ctype.tbl.norm = meas.ctype.tbl / rowSums(meas.ctype.tbl)

  # export num cell table
  num.cell.df = as.data.frame(samp.ctype.tbl)
  colnames(num.cell.df) = paste0('num_', col.to.ctype[colnames(num.cell.df)])
  num.cell.df.fil = cbind(data.frame(indiv_id=names(indiv.to.used.samp)), num.cell.df[indiv.to.used.samp,])
  write.csv(num.cell.df.fil, file=file.path(SUPP.TABLE.DIR, 'num_cells_per_ctype.csv'), quote=F, row.names=F)

  # edf 7B
  replicate.fig.dir = file.path(fig.dir, 'replicate_cell_types')
  dir.create(replicate.fig.dir, showWarnings=F)
  for (i in 1:ncol(samp.ctype.tbl.norm)) {
    cur.col = colnames(samp.ctype.tbl.norm)[i]
    xvals1 = meas.ctype.tbl.norm[unlist(tech.reps[should.use.tech.rep,1]), cur.col]
    yvals1 = meas.ctype.tbl.norm[unlist(tech.reps[should.use.tech.rep,2]), cur.col]
    xvals2 = samp.ctype.tbl.norm[unlist(bio.reps[should.use.bio.rep,1]), cur.col]
    yvals2 = samp.ctype.tbl.norm[unlist(bio.reps[should.use.bio.rep,2]), cur.col]
    cur.lim = range(c(xvals1, yvals1, xvals2, yvals2))
    png(file.path(replicate.fig.dir, paste0(cur.col, '.png')), width=1100, height=500)
    par(mfrow=c(1, 2))
    plot(xvals1, yvals1, pch=21, bg=cur.col, main=col.to.ctype[cur.col], cex=2, xlim=cur.lim, ylim=cur.lim)
    abline(b=1, a=0, col=2, lwd=2)
    plot(xvals2, yvals2, pch=21, bg=cur.col, main=col.to.ctype[cur.col], cex=2, xlim=cur.lim, ylim=cur.lim)
    abline(b=1, a=0, col=2, lwd=2)
    dev.off()
  }

  # main replicate figures
  # fig 2C
  for (i in 1:3) {
    cur.name = c('clp', 'mebempe', 'mebempl')[i]
    cur.types = list(c('blue', 'darkblue', '#6E1C71', 'lightblue'), c('brown', 'gold', '#7F9D00'), c('purple', 'red', 'plum', '#EEBB6E'))[[i]]
    xvals = rowSums(samp.ctype.tbl.norm[unlist(bio.reps[should.use.bio.rep,1]), cur.types])
    yvals = rowSums(samp.ctype.tbl.norm[unlist(bio.reps[should.use.bio.rep,2]), cur.types])
    cur.lim = range(c(xvals, yvals))
    png(file.path(replicate.fig.dir, paste0(cur.name, '.png')), width=500, height=500)
    plot(xvals, yvals, pch=21, bg=1, main=cur.name, cex=2, xlim=cur.lim, ylim=cur.lim)
    abline(b=1, a=0, col=2, lwd=2)
    dev.off()
  }

  # clinical correlations
  # fig 2E, fig 3A, fig 3B, fig 3C, edf 8A
  clinical.vs.cell.type.freq(indiv.to.used.samp.fil, samp.ctype.tbl, ages, sexes)

  # N122
  # edf 8B
  analyze.n122(our.assign, our.annotation, our.legc, fig.dir=fig.dir)

  # rdw vs clonality in a separate cohort
  # fig 2F
  fig2.rdw.vs.mut(fig.dir)

  # analyzing CD34+ cell frequency in another dataset
  # fig 3D
  garvan.analysis()

  our.cdata.x = our.cdata$X
  # fig 3E, fig 3F, fig 3G, fig 3H
  # edf 8C, edf 8D, edf 8E, edf 8F , edf 8G
  de.analysis(our.cdata, our.cdata.x, our.mdata, our.legc, our.annotation, our.assign, cell.to.samp, indiv.to.used.samp, indiv.to.used.samp, fig.dir=fig.dir)

  # fig 3I, fig 3J, fig 3K, edf 8H
  sync.analysis(our.annotation, our.assign, our.cdata, our.legc, our.cdata.x, indiv.to.used.samp, cell.to.meas, cell.to.samp, samp.ctype.tbl)

  create.supp.tables(our.cdata, our.annotation, our.assign, indiv.to.used.samp)

}

plot.main.traj <- function(our.legc, our.annotation, cell.to.samp, indiv.to.used.samp.fil, ages, indivs.cols) {
  # big frequency figure
  mebemp.cells = names(our.assign)[our.assign %in% names(our.annotation)[our.annotation %in% c('gold', '#EEBB6E', 'plum')]]
  mebemp.avp.scores = our.legc['AVP', our.assign[mebemp.cells]]
  names(mebemp.avp.scores) = mebemp.cells
  mebemp.bin = sort(rep(1:10, length.out=length(mebemp.avp.scores)))
  names(mebemp.bin) = names(mebemp.avp.scores)[order(mebemp.avp.scores, decreasing=T)]
#
  clp.cells = names(our.assign)[our.assign %in% names(our.annotation)[our.annotation %in% c('blue', 'darkblue')]]
  clp.avp.scores = our.legc['AVP', our.assign[clp.cells]]
  names(clp.avp.scores) = clp.cells
  clp.bin = sort(rep(2:4, length.out=length(clp.avp.scores)))
  names(clp.bin) = names(clp.avp.scores)[order(clp.avp.scores, decreasing=T)]
  clp.bin = -clp.bin

  hsc.cells = names(our.assign)[our.assign %in% names(our.annotation)[our.annotation == 'brown']]
  clpe.cells = names(our.assign)[our.assign %in% names(our.annotation)[our.annotation == '#6E1C71']]
  hsc.bin = rep(0, length(hsc.cells))
  names(hsc.bin) = hsc.cells
  clpe.bin = rep(-1, length(clpe.cells))
  names(clpe.bin) = clpe.cells
  main.bin = c(mebemp.bin, hsc.bin, clpe.bin, clp.bin)

  freq.tbl = table(main.bin, cell.to.samp[names(main.bin)])
  samps.to.plot = intersect(colnames(freq.tbl)[colSums(freq.tbl) > 500], indiv.to.used.samp.fil)
  freq.tbl.norm = t(t(freq.tbl[,samps.to.plot]) / colSums(freq.tbl[,samps.to.plot]))
  #freq.tbl.legc = log2(1e-3 + freq.tbl.norm)
  freq.tbl.relative = log2(1e-3 + freq.tbl.norm / apply(freq.tbl.norm, 1, median))
  cur.samps.ord.orig = names(sort(colSums(freq.tbl.relative * 1:nrow(freq.tbl.relative))))
  tmp.early.scores = colSums(freq.tbl[as.character(1:5),]) / colSums(freq.tbl[as.character(1:10),])
  tmp.clp.scores = colMeans(freq.tbl.relative[1:3,])
  clp.splitted = split(tmp.clp.scores, cut(tmp.clp.scores, c(-1000, -0.5, 0.5, 1000)))

  cur.samps.ord = c(names(clp.splitted[[1]])[order(tmp.early.scores[names(clp.splitted[[1]])], decreasing=T)], 
                     names(clp.splitted[[2]])[order(tmp.early.scores[names(clp.splitted[[2]])], decreasing=T)], 
                     names(clp.splitted[[3]])[order(tmp.early.scores[names(clp.splitted[[3]])], decreasing=T)])

  samp.freq.cls = rep(NA, length(samps.to.plot))
  names(samp.freq.cls) = samps.to.plot
  samp.freq.cls[samps.to.plot %in% names(clp.splitted[[1]])] = 1
  samp.freq.cls[samps.to.plot %in% names(clp.splitted[[2]])] = 2
  samp.freq.cls[samps.to.plot %in% names(clp.splitted[[3]])] = 3
  samp.freq.cls = samp.freq.cls * 2
  samp.freq.cls = samp.freq.cls - ifelse(tmp.early.scores[samps.to.plot] > 0.5, 1, 0)
  samp.freq.cls = 7 - samp.freq.cls

  #cur.samps.ord = cur.samps.ord[cur.samps.ord %in% samps.to.plot]
  cur.indivs.ord = sapply(strsplit(cur.samps.ord, '-'), function(x) x[2])
  tmp.df = data.frame(samps_cols=indivs.cols[cur.indivs.ord], age=pmin(pmax(ages[cur.indivs.ord], 20), 85), 
                      early_score=tmp.early.scores[cur.samps.ord], freq_cls=samp.freq.cls[cur.samps.ord])
  rownames(tmp.df) = cur.samps.ord
  tmp.cols = unique(indivs.cols)
  names(tmp.cols) = tmp.cols
  age.cols = colorRampPalette(c('white', 'grey', 'black'))(66)
  names(age.cols) = 20:85
  #cls.cols = RColorBrewer::brewer.pal(7, 'Set3')[2:7]
  cls.cols = rev(RColorBrewer::brewer.pal(7, 'Set3')[c(1, 3:7)])
  names(cls.cols) = 1:6

  #pheatmap(t(pmax(pmin(freq.tbl.relative[,cur.samps.ord], 3), -3))[,15:1], breaks=seq(-3, 3, length.out=length(shades)), 
  #         col=shades, cluster_cols=F, cluster_rows=F, annotation_row=tmp.df, annotation_colors=list(samps_cols=tmp.cols, age=age.cols))
  shades = colorRampPalette(c("darkblue", "blue", "white",'red', "darkred"))(101)
  # fig 2D bottom
  pheatmap(t(pmax(pmin(freq.tbl.relative[,cur.samps.ord], 1), -1))[,15:1], breaks=seq(-1, 1, length.out=length(shades)), 
           col=shades, cluster_cols=F, cluster_rows=F, annotation_row=tmp.df, annotation_colors=list(samps_cols=tmp.cols, age=age.cols, freq_cls=cls.cols),
	   height=14, width=7.5, filename=file.path(fig.dir, 'freq_heatmap.png'))

  # specific examples
  bg.samps = samps.to.plot
  indivs.to.plot.alone = c('N115', 'N165', 'N84', 'N261', 'N107', 'N179')
  samps.to.plot.alone = indiv.to.used.samp.fil[indivs.to.plot.alone]
  median.freqs = apply(freq.tbl.norm[,bg.samps], 1, median)
  q95.freqs = apply(freq.tbl.norm[,bg.samps], 1, quantile, 0.95)
  q05.freqs = apply(freq.tbl.norm[,bg.samps], 1, quantile, 0.05)
  #par(mfrow=c(3, 2), mar=c(2, 2, 2, 2))
  #png(file.path(fig.dir, 'selected_indivs_freq.png'), width=3000, height=700)
  png(file.path(fig.dir, 'selected_indivs_freq.png'), width=7000, height=600)
  #par(mfrow=c(2, 3), mar=c(6, 2, 2, 2))
  par(mfrow=c(1, 6), mar=c(6, 2, 2, 2))
  # fig 2D top
  for (i in seq_along(samps.to.plot.alone)) {
    plot(rev(as.numeric(rownames(freq.tbl.norm))), median.freqs, ylim=c(0, 0.2), type='l', lwd=8, main=indivs.to.plot.alone[i], lty='dashed', yaxt='n', xaxt='n', xlab='')
    cur.labels = rep('', 15)
    #cur.labels = c(paste0('MEBEMP-', 10:1), 'HSC', paste0('CLP-', 1:4))
    #axis(side=1, at=-4:10, labels=cur.labels, las=2)
    lines(rev(as.numeric(rownames(freq.tbl.norm))), q95.freqs, ylim=c(0, 0.2), type='l', lwd=14, main=indivs.to.plot.alone[i], lty='dashed', col='grey')
    lines(rev(as.numeric(rownames(freq.tbl.norm))), q05.freqs, ylim=c(0, 0.2), type='l', lwd=14, main=indivs.to.plot.alone[i], lty='dashed', col='grey')
    lines(as.numeric(rownames(freq.tbl.norm)), freq.tbl.norm[nrow(freq.tbl.norm):1,samps.to.plot.alone[i]], 
          lwd=35, ylim=c(0, 0.2), col=cls.cols[samp.freq.cls[samps.to.plot.alone[i]]])
    box(lwd=3)
  }
  dev.off()
}

get.ctype.distrib <- function(our.assign, our.annotation, cell.to.val, used.mcs=NULL) {
  mc_pat = t(table(cell.to.val[names(our.assign)], our.assign))
  mc_pat = mc_pat[-1,]
  mpat_f = matrix(mc_pat, ncol=ncol(mc_pat), dimnames=dimnames(mc_pat))
  stopifnot(all(rownames(mpat_f) %in% names(our.annotation)))
  if (is.null(used.mcs)) {
    pat_type = t(tgstat::tgs_matrix_tapply(t(mpat_f), our.annotation[rownames(mpat_f)], sum))
  } else {
    pat_type = t(tgstat::tgs_matrix_tapply(t(mpat_f[used.mcs,]), our.annotation[used.mcs], sum))
  }
  return(pat_type)
}

plot.cell.type.freqs <- function(our.assign, our.annotation, our.legc, cell.to.samp, indiv.to.used.samp.fil, fig.dir=file.path(BASE.FIG.DIR, 'fig2')) {
  
  cd34.pos.mcs = colnames(our.legc)[our.legc['CD34',] > -14.5]
  pat.type = get.ctype.distrib(our.assign, our.annotation, cell.to.samp, cd34.pos.mcs)
  
  cur.used.ctypes = c('purple', 'red', 'plum', '#EEBB6E', "#7F9D00", 'gold', 'brown', "#6E1C71", 'blue', 'darkblue', 'lightblue')
  pat.type.fil = pat.type[,cur.used.ctypes]
  high.cov.samps = rownames(pat.type.fil)[rowSums(pat.type.fil) >= 1000]
  used.samps = intersect(high.cov.samps, indiv.to.used.samp.fil)
  pat_type_n1 = pat.type.fil[used.samps,] / rowSums(pat.type.fil[used.samps,])
  
  pheatmap::pheatmap(tgs_cor(pat_type_n1), cluster_rows=F, cluster_cols=F, color=colorRampPalette(c("darkblue","white","darkred"))(1000), breaks=c(-1,seq(-0.7,0.7,l=999),1),file=file.path(fig.dir, "pat_type_cor.png"), w=8,h=8)
  
  pheatmap::pheatmap(tgs_cor(pat_type_n1, spearman=T), cluster_rows=F, cluster_cols=F, color=colorRampPalette(c("darkblue","white","darkred"))(1000), breaks=c(-1,seq(-0.7,0.7,l=999),1),file=file.path(fig.dir, "pat_type_scor.png"), w=8,h=8)
  
  # fig 2B
  png(file.path(fig.dir, "pat_type_box.png"), w=700, h=400)
  boxplot(log10(pat_type_n1+1e-2),las=2,col=colnames(pat_type_n1), boxwex=0.8, lwd=2)
  dev.off()

  print(round(colMeans(pat_type_n1), 3))
  print(round(apply(pat_type_n1, 2, sd), 3))

}

de.between.indivs.analysis <- function(our.cdata, our.mdata, our.annotation, our.assign, indiv.to.used.samp.fil, cell.to.samp, fig.dir=BASE.FIG.DIR) {
  # separately for CLPs and MEBEMPs
  cd34.mebemp.ctypes = c("purple", "red", "plum", "#EEBB6E", "#7F9D00", "gold")
  cd34.clp.ctypes = c("#6E1C71", "blue", "darkblue", "lightblue")
  sex.genes = c(read.table(SEX.GENES.PATH, stringsAsFactors=F)[,1], 'NOVA1', 'XIST')
  our.cdata.x = our.cdata$X
  our.mdata.x = our.mdata$layers[['total_umis']]
  #rownames(our.mdata.x) = paste0('mc', 0:(nrow(our.mdata.x) - 1))
  rownames(our.mdata.x) = paste0('mc', 0:(nrow(our.mdata.x) - 1))
  our.annotation.ext = our.annotation
  our.annotation.ext['mc-1'] = 'outlier'

  for (i in 1:2) {
    cur.ctypes = list(cd34.mebemp.ctypes, cd34.clp.ctypes)[[i]]
    cur.cells = names(our.assign)[our.annotation.ext[as.character(our.assign)] %in% cur.ctypes]

    set.seed(42)

    cur.cells.fil = c()
    cur.exp.pooled.list = list()
    tmp.samps = unique(cell.to.samp[cur.cells])
    for (j in 1:length(tmp.samps)) {
      cur.samp = tmp.samps[j]
      cells = cur.cells[cell.to.samp[cur.cells] == cur.samp]
      print(unique(cell.to.samp[cells]))
      print(Sys.time())
      cur.samp.mat = our.cdata.x[cells,,drop=F]
      cur.mat.ds = scm_downsamp(t(cur.samp.mat), 500)
      cur.exp.pooled.list[[j]] = rowSums(cur.mat.ds)
      cur.cells.fil = c(cur.cells.fil, colnames(cur.mat.ds))
    }
    cur.exp.pooled = do.call(cbind, cur.exp.pooled.list)
    colnames(cur.exp.pooled) = tmp.samps

    cur.thresh = c(3e5, 1e5)[i]
    cur.samps = colnames(cur.exp.pooled)[colSums(cur.exp.pooled) > cur.thresh]
    cur.exp.pooled.frac = t(t(cur.exp.pooled[,cur.samps]) / colSums(cur.exp.pooled[,cur.samps]))
    cur.indiv.legc = log2(1e-5 + cur.exp.pooled.frac)

    cur.kruskal.pvals.df = read.table(BATCH.EFFECT.PVALS.PATH, header=T)
    cur.kruskal.pvals = cur.kruskal.pvals.df[,2]
    names(cur.kruskal.pvals) = cur.kruskal.pvals.df[,1]
    cur.non.batchy.genes = setdiff(rownames(cur.exp.pooled.frac), names(cur.kruskal.pvals)[cur.kruskal.pvals < 1e-4])
  
    tmp.tbl = table(our.assign[cur.cells.fil])
    our.mdata.x.ds = t(scm_downsamp(t(our.mdata.x), 9e4))
    cur.exp.pooled.expected = do.call(cbind, tapply(cur.cells.fil, cell.to.samp[cur.cells.fil], function(cells) {colSums(our.mdata.x.ds[our.assign[cells],,drop=F])}))

    cur.exp.pooled.expected.frac = t(t(cur.exp.pooled.expected[,cur.samps]) / colSums(cur.exp.pooled.expected[,cur.samps]))
    cur.indiv.expected.legc = log2(1e-5 + cur.exp.pooled.expected.frac)

    # remove one very batchy experiment
    non.batchy.samps = intersect(setdiff(cur.samps, grep('demux_07_02_22-', cur.samps, v=T)), indiv.to.used.samp.fil)
    tmp.obs = cur.indiv.legc[,non.batchy.samps]
    tmp.exp = cur.indiv.expected.legc[,non.batchy.samps]
    tmp.diff = (tmp.obs - tmp.exp)[cur.non.batchy.genes,]

    min.expression = c(-14.5, -13.5)[i]
    is.expressed = (pmax(tmp.obs, tmp.exp) > min.expression)[cur.non.batchy.genes,]
    diff.genes = rownames(tmp.diff)[rowSums(abs(tmp.diff) > 1 & is.expressed) > 1]
    diff.genes.no.sex = setdiff(diff.genes, sex.genes)
    if (i == 1) {
      # removing JUN, histones, CLP genes, and genes that are batch specific (RNF11-PRKAR2B, CHMP1A-BTBD2)
      # removing also B module (IGHM-LTB), histones, and CHMP1A-BTBD2 is another module
      # The genes that are seen in the heatmap are the interferon response
      diff.genes.no.sex = setdiff(diff.genes.no.sex, c(c('JUN', 'FOSB', 'DUSP1', 'FOS', 'RNF11', 'C2orf88', 'RUFY1', 'ARHGAP18', 'PRKAR2B'),
                                  grep('HIST', diff.genes.no.sex, v=T),
				  c('IGHM', 'CD79A', 'IGHA1', 'JCHAIN', 'LTB', 'CHMP1A', 'PTOV1', 'TARS2', 'ESRRA', 'INTS1', 'EHD2', 'AP3M1', 'TESPA1', 'PREX1', 'BTBD2')))
    } else {
      # revmoving jun etc
      diff.genes.no.sex = setdiff(diff.genes.no.sex, c('JUN', 'FOSB', 'DUSP1', 'FOS'))
    }
    shades = colorRampPalette(c("darkblue", "blue", "white","red", "darkred"))(100)
    shades[40:60] = 'white'
    # edf 4A, 4B
    pheatmap(pmax(pmin(tmp.diff[diff.genes.no.sex,], 3), -3), col=shades, breaks=seq(-3, 3, length.out=length(shades) + 1), filename=file.path(fig.dir, sprintf('de_screen_%s.png', i)), height=10, width=5)
  }
}

get.biological.replicates <- function(all.samp, indiv.to.used.samp) {
  samp.to.indiv = sapply(strsplit(all.samp, '-'), function(x) x[2])
  names(samp.to.indiv) = all.samp
  tmp.tbl = table(samp.to.indiv)
  indivs.with.replicates = names(tmp.tbl)[tmp.tbl > 1]
  all.rep1 = c()
  all.rep2 = c()
  for (cur.indiv in indivs.with.replicates) {
    cur.samples = names(samp.to.indiv)[samp.to.indiv == cur.indiv]
    cur.samples.for = unique(gsub('_ultima', '', cur.samples))
    if (length(cur.samples.for) == 1) {
      next
    }
    is.used.samp = cur.samples %in% indiv.to.used.samp
    stopifnot(sum(is.used.samp) == 1)

    used.samp = cur.samples[is.used.samp] 
    if (!(used.samp %in% cur.samples.for)) {
      tmp.sample.pref = sapply(strsplit(cur.samples.for, '_|-'), function(x) paste(x[1:4], collapse='_'))
      used.sample.pref = sapply(strsplit(used.samp, '_|-'), function(x) paste(x[1:4], collapse='_'))
      cur.samples.for[used.sample.pref == tmp.sample.pref] = used.samp
    }
    is.used.samp.for = cur.samples.for %in% indiv.to.used.samp
    stopifnot(sum(is.used.samp.for) == 1)
    stopifnot(length(is.used.samp.for) == 2)
    samp.prefix = sapply(strsplit(cur.samples.for, '-'), function(x) x[1])
    tmp.ord = order.samp(samp.prefix)
    if (any(is.used.samp.for[tmp.ord] != c(T, F))) {
     print(sprintf('Not using individual %s because good sample is not the first sample', cur.indiv))
     next
    }
    stopifnot(all(cur.samples.for[tmp.ord] %in% all.samp))

    print('**************')
    print(cur.indiv)
    print(samp.prefix[tmp.ord])
    print('**************')
    all.rep1 = c(all.rep1, cur.samples.for[tmp.ord][1])
    all.rep2 = c(all.rep2, cur.samples.for[tmp.ord][2])
  }
  bio.reps = data.frame(all.rep1, all.rep2)
}

order.samp <- function(meas) {
  meas.no.indiv = sapply(strsplit(meas, '-'), function(x) x[1])
  tmp.splitted = strsplit(meas.no.indiv, '_')
  tmp.year.char = sapply(tmp.splitted, function(x) x[4])
  tmp.month.char = sapply(tmp.splitted, function(x) x[3])
  tmp.day.char = sapply(tmp.splitted, function(x) x[2])
  stopifnot(all(nchar(tmp.year.char) == 2))
  stopifnot(all(nchar(tmp.month.char) == 2))
  stopifnot(all(nchar(tmp.day.char) == 2))
  tmp.year = as.numeric(tmp.year.char)
  tmp.month = as.numeric(tmp.month.char)
  tmp.day = as.numeric(tmp.day.char)
  stopifnot(all(!is.na(tmp.year)))
  stopifnot(all(!is.na(tmp.month)))
  stopifnot(all(!is.na(tmp.day)))
  return(order(tmp.year, tmp.month, tmp.day))
}

get.indiv.to.vaf <- function(indiv.to.used.samp, indiv.id.map=NULL) {
  if (is.null(indiv.id.map)) {
    indiv.id.map = get.indiv.id.map()
  }
  all.muts.df = get.mutations.new(indiv.to.used.samp, indiv.id.map)
  indiv.to.vaf = rep(0, length(indiv.to.used.samp))
  names(indiv.to.vaf) = names(indiv.to.used.samp)
  max.vafs = tapply(all.muts.df$Avg_VAF, all.muts.df$sample, max)
  max.vafs.fil = max.vafs[names(max.vafs) %in% all.indivs]
  indiv.to.vaf[names(max.vafs.fil)] = max.vafs.fil
  return(indiv.to.vaf)
}

permut_test = function(pats, feat, label, pat_type_n1, win=2, num.perm=1000) {
  rstats = c()
  cs = apply(pat_type_n1[pats,], 2, cor, feat, use="pairwise.complete.obs", m="spearman")
  cs2 = rollmean(cs,win)
  stat = max(abs(cs2))
  
  set.seed(42)
  for(i in 1:num.perm) {
  	rpats = sample(pats,length(pats))
  	rcs = apply(pat_type_n1[rpats,], 2, cor, feat, use="pairwise.complete.obs", m="spearman")
  	rcs2 = rollmean(rcs,win)
  	rstats = c(rstats, max(abs(rcs2)))
  }
  #png(sprintf("figs/pat_typecbc/%s.png", label), w=800,h=400)
  ylim = c(-max(abs(cs)), max(abs(cs)))
  barplot(cs, col=colnames(pat_type_n1), las=2, main=sprintf("%s, N=%s, %d/%d",label, length(pats), sum(rstats>=stat), num.perm), ylim=ylim)
  return(list(sum(rstats >= stat), num.perm, stat, rstats))
}

clinical.vs.cell.type.freq <- function(indiv.to.used.samp.fil, samp.ctype.tbl, ages, sexes, fig.dir=file.path(BASE.FIG.DIR, 'fig2')) {
  indiv.id.map = get.indiv.id.map()
  cbc.new = get.clinical.values.new(indiv.id.map)
  stopifnot(all(sort(names(indiv.to.used.samp.fil)[which(!(names(indiv.to.used.samp.fil) %in% rownames(cbc.new)))]) == c('N170', 'N307')))

  used.samp.to.indiv = names(indiv.to.used.samp.fil)
  names(used.samp.to.indiv) = indiv.to.used.samp.fil

  samp.ctype.tbl.norm = samp.ctype.tbl / rowSums(samp.ctype.tbl)
  high.cov.samps = intersect(rownames(samp.ctype.tbl)[rowSums(samp.ctype.tbl) >= 500], indiv.to.used.samp.fil)
  high.cov.indivs = used.samp.to.indiv[high.cov.samps]

  indivs.for.cbc = setdiff(high.cov.indivs, c('N170', 'N307'))
  samps.for.cbc = indiv.to.used.samp.fil[indivs.for.cbc]

  pat_type_n1 = samp.ctype.tbl.norm

  # all samples
  # fig 2E
  #par(mfrow=c(4, 5))
  cur.dir = file.path(fig.dir, 'cbc_all')
  dir.create(cur.dir, showWarnings=F)
  for (i in 1:20) {
    png(file.path(cur.dir, paste0(i, '.png')), height=500, width=800)
    permut_test(indiv.to.used.samp.fil[indivs.for.cbc], cbc.new[indivs.for.cbc, i], colnames(cbc.new)[i], pat_type_n1, win=3)
    dev.off()
  }

  # males
  #par(mfrow=c(4, 5))
  indiv.to.used.samp.male = indiv.to.used.samp.fil[indivs.for.cbc][sexes[indivs.for.cbc] == 'male']
  cur.dir = file.path(fig.dir, 'cbc_male')
  dir.create(cur.dir, showWarnings=F)
  for (i in 1:20) {
    png(file.path(cur.dir, paste0(i, '.png')), height=500, width=800)
    permut_test(indiv.to.used.samp.male, cbc.new[names(indiv.to.used.samp.male), i], colnames(cbc.new)[i], pat_type_n1, win=3)
    dev.off()
  }

  # females
  #par(mfrow=c(4, 5))
  indiv.to.used.samp.female = indiv.to.used.samp.fil[indivs.for.cbc][sexes[indivs.for.cbc] == 'female']
  cur.dir = file.path(fig.dir, 'cbc_female')
  dir.create(cur.dir, showWarnings=F)
  for (i in 1:20) {
    png(file.path(cur.dir, paste0(i, '.png')), height=500, width=800)
    permut_test(indiv.to.used.samp.female, cbc.new[names(indiv.to.used.samp.female), i], colnames(cbc.new)[i], pat_type_n1, win=3)
    dev.off()
  }


  # mutation analysis
  indiv.to.vaf = get.indiv.to.vaf(indiv.to.used.samp.fil, indiv.id.map)
  excluded.indivs = c('N116', 'N130', 'N180', 'N181', 'N201')
  excluded.samps = indiv.to.used.samp.fil[names(indiv.to.used.samp.fil) %in% excluded.indivs]
  sel.samps.unfil = intersect(rownames(samp.ctype.tbl)[rowSums(samp.ctype.tbl) > 200], indiv.to.used.samp.fil)
  sel.samps = setdiff(sel.samps.unfil, excluded.samps)

  clp.fracs = rowSums(samp.ctype.tbl.norm[sel.samps, c('blue', 'darkblue', 'lightblue', '#6E1C71')])
  clp.fracs.splitted = split(clp.fracs, indiv.to.vaf[used.samp.to.indiv[sel.samps]] > 0)
  print(wilcox.test(clp.fracs.splitted[[1]], clp.fracs.splitted[[2]]))
  # edf 8A
  png(file.path(fig.dir, 'clp_mut_boxplot.png'), width=300, height=600)
  boxplot(clp.fracs.splitted, lwd=3)
  dev.off()

  mut.sel.samps = sel.samps[indiv.to.vaf[used.samp.to.indiv[sel.samps]] > 0]
  png(file.path(fig.dir, 'clp_mut_scatter.png'), width=600, height=600)
  plot(clp.fracs[mut.sel.samps], indiv.to.vaf[used.samp.to.indiv[mut.sel.samps]], pch=21, bg=indivs.cols[used.samp.to.indiv[mut.sel.samps]], cex=3)
  dev.off()


  # age correlations
  used.samps = setdiff(high.cov.samps, high.cov.samps[indiv.to.vaf[high.cov.indivs] > 0])

  # age x sex boxplots
  tmp.young.males = used.samps[ages[samp.to.indiv[used.samps]] <= 50 & sexes[samp.to.indiv[used.samps]] == 'male']
  tmp.old.males = used.samps[ages[samp.to.indiv[used.samps]] > 60 & sexes[samp.to.indiv[used.samps]] == 'male']
  tmp.young.females = used.samps[ages[samp.to.indiv[used.samps]] <= 50 & sexes[samp.to.indiv[used.samps]] == 'female']
  tmp.old.females = used.samps[ages[samp.to.indiv[used.samps]] > 60 & sexes[samp.to.indiv[used.samps]] == 'female']

  #tmp.cols = c('lightgreen', '#FFD580', RColorBrewer::brewer.pal(5, 'Set1')[3], RColorBrewer::brewer.pal(5, 'Set1')[5])
  tmp.cols = c('blue', 'red', 'blue', 'red')
  clp.freqs = rowSums(samp.ctype.tbl.norm[used.samps, c('blue', 'darkblue', '#6E1C71')])
  mebemp.freqs = rowSums(samp.ctype.tbl.norm[used.samps, c('#EEBB6E', 'plum', 'red', 'purple')])
  cur.clp.list = list(clp.freqs[tmp.young.males], clp.freqs[tmp.young.females], clp.freqs[tmp.old.males], clp.freqs[tmp.old.females])
  cur.mebemp.list = list(mebemp.freqs[tmp.young.males], mebemp.freqs[tmp.young.females], mebemp.freqs[tmp.old.males], mebemp.freqs[tmp.old.females])
  print('kruskal clp')
  print(kruskal.test(cur.clp.list))
  print('kruskal mebemp')
  print(kruskal.test(cur.mebemp.list))
  # fig 3A
  png(file.path(fig.dir, 'age_sex_clp_box.png'), height=800, width=400)
  boxplot(cur.clp.list, col=tmp.cols, lwd=3)
  dev.off()
  png(file.path(fig.dir, 'age_sex_mebemp_box.png'), height=800, width=400)
  boxplot(cur.mebemp.list, col=tmp.cols, lwd=3)
  dev.off()

  mebempl.mpp.ratio = log2(samp.ctype.tbl.norm[used.samps, 'plum'] / samp.ctype.tbl.norm[used.samps, 'gold'])
  mebempl.mpp.ratio.list = list(mebempl.mpp.ratio[tmp.young.males], mebempl.mpp.ratio[tmp.young.females], mebempl.mpp.ratio[tmp.old.males], mebempl.mpp.ratio[tmp.old.females])
  print('kruskal mebempl mpp ratio')
  print(kruskal.test(mebempl.mpp.ratio.list))
  # fig 3B
  png(file.path(fig.dir, 'age_sex_mebempl_mpp_ratio_box.png'), height=800, width=400)
  boxplot(mebempl.mpp.ratio.list, col=tmp.cols, lwd=3)
  dev.off()

  hsc.freqs = samp.ctype.tbl.norm[used.samps, 'brown']
  cur.hsc.list = list(hsc.freqs[tmp.young.males], hsc.freqs[tmp.young.females], hsc.freqs[tmp.old.males], hsc.freqs[tmp.old.females])
  print('kruskal hsc')
  print(kruskal.test(cur.hsc.list))
  # fig 3C
  png(file.path(fig.dir, 'age_sex_hsc_box.png'), height=800, width=400)
  boxplot(cur.hsc.list, col=tmp.cols, lwd=3)
  dev.off()

}

get.mutations.new <- function(indiv.to.used.samp, indiv.id.map) {
  arch3.muts.path = file.path(MUTATION.FILES.DIR, 'arch3.csv')
  arch4.muts.path = file.path(MUTATION.FILES.DIR, 'arch4.csv')
  arch3.muts = read.csv(arch3.muts.path, stringsAsFactors=F)
  arch4.muts = read.csv(arch4.muts.path, stringsAsFactors=F)

  sample.id.regex = regexpr('\\d+', arch4.muts$Sample_num)
  tmp.nums = as.numeric(substr(arch4.muts$Sample_num, sample.id.regex, sample.id.regex + attr(sample.id.regex, 'match.length') - 1))
  old.indivs.in.arch4 = unique(tmp.nums[tmp.nums < 185])
  stopifnot(sum(paste0('N', old.indivs.in.arch4) %in% names(indiv.to.used.samp)) == 0)

  trans.func = function(cur.sample) {
    if (cur.sample %in% names(indiv.id.map)) {
      return(indiv.id.map[cur.sample])
    } else {
      return(cur.sample)
    }
  }

  arch3.muts$sample = sapply(paste0('N', arch3.muts$Sample_num), trans.func)
  arch4.muts$sample = sapply(arch4.muts$Sample_num, trans.func)
  stopifnot(length(Reduce(intersect, list(arch3.muts$sample, arch4.muts$sample, names(indiv.to.used.samp)))) == 0)

  is.na.avg.vaf = is.na(arch3.muts$Avg_VAF)
  arch3.muts$Avg_VAF[is.na.avg.vaf] = arch3.muts$VAF[is.na.avg.vaf]

  # assert that an individual doesn't appear twice with different ids in the same file
  ids.per.indiv = colSums(table(arch4.muts$Sample_num, arch4.muts$sample) > 0)
  indivs.with.many.ids = intersect(names(ids.per.indiv)[ids.per.indiv > 1], names(indiv.to.used.samp))
  stopifnot(length(indivs.with.many.ids) == 0)

  rel.cols = c('Gene.refGene', 'Chr', 'Pos', 'Ref', 'Alt', 'Avg_VAF', 'sample')
  stopifnot(all(rel.cols %in% colnames(arch3.muts)))
  stopifnot(all(rel.cols %in% colnames(arch4.muts)))
  all.muts.df = rbind(arch3.muts[,rel.cols], arch4.muts[,rel.cols])
  stopifnot(all(!is.na(all.muts.df$Avg_VAF)))
  return(all.muts.df)
}

get.clinical.values.from.file <- function(clin.fname) {
  long.cbc = read.csv(clin.fname, stringsAsFactors=F)[,1:22]
  cur.num = long.cbc[1, 1]
  for (i in 1:nrow(long.cbc)) {
    if (is.na(long.cbc[i, 1])) {
      long.cbc[i, 1] = cur.num
    } else if (long.cbc[i, 1] != cur.num) {
        cur.num = long.cbc[i, 1]
    }
    long.cbc[i, 1] = cur.num
  }
  is.empty.row = rowSums(is.na(long.cbc[,3:22])) == 20
  long.cbc = long.cbc[!is.empty.row,]
  stopifnot(all(grepl('\\d\\d\\.\\d\\d\\.\\d\\d', long.cbc$Date)))
  long.cbc$day = substr(long.cbc$Date, 1, 2)
  long.cbc$month = substr(long.cbc$Date, 4, 5)
  long.cbc$year = substr(long.cbc$Date, 7, 8)
  long.cbc$full_year = ifelse(as.numeric(long.cbc$year) > 25, paste0('19', long.cbc$year), paste0('20', long.cbc$year))
  long.cbc$date_str = sprintf('%s/%s/%s', long.cbc$day, long.cbc$month, long.cbc$full_year)
  dates = as.Date(long.cbc$date_str, format='%d/%m/%Y')
  tmp.rle = rle(long.cbc[,1])
  cur.start.index = 1
  for (i in 1:length(tmp.rle$lengths)) {
    cur.end.index = sum(tmp.rle$lengths[1:i])
    if (!(all(sort(dates[cur.start.index:cur.end.index], decreasing=T) == dates[cur.start.index:cur.end.index]))) {
      print(tmp.rle$values[i])
      #print(i)
      #print(sort(dates[cur.start.index:cur.end.index], decreasing=T))
      #print(dates[cur.start.index:cur.end.index])
    }
    cur.start.index = cur.start.index + tmp.rle$length[i]
  }

  # translate to quantiles and make sure they are retained
  colnames(long.cbc)[3:22] = c('WBC', 'RBC', 'Hemoglobin', 'Hematocrit', 'MCV', 'MCH', 'MCHC', 'RDW', 'Platelets', 'MPV', 
                    'Neutro%', 'Lympho%', 'Mono%', 'Eosinophils%', 'Basophils%', 'Neutro#', 'Lympho#', 'Mono#', 'Eosinophils#', 'Basophils#')

  long.cbc.recent = filter(long.cbc, full_year >= 2018)
  indivs.with.cbc = as.character(unique(long.cbc.recent[,1]))
  corrected.cbc = do.call(rbind, lapply(indivs.with.cbc, function(cur.indiv) {
    #colMeans(long.cbc.recent[long.cbc.recent[,1] == cur.indiv, 3:22, drop=F], na.rm=T)
    apply(long.cbc.recent[long.cbc.recent[,1] == cur.indiv, 3:22, drop=F], 2, median, na.rm=T)
  }))
  rownames(corrected.cbc) = indivs.with.cbc
  return(corrected.cbc)
}

get.clinical.values.new <- function(indiv.id.map) {

  clinical.values.prev = get.clinical.values.from.file(PREV.CLIN.VALUES.PATH)
  clinical.values.new = get.clinical.values.from.file(NEW.CLIN.VALUES.PATH)
  rownames(clinical.values.prev) = paste0('N', rownames(clinical.values.prev))
  rownames(clinical.values.new) = paste0('N', rownames(clinical.values.new))
  stopifnot(length(intersect(c(rownames(clinical.values.new), indiv.id.map[rownames(clinical.values.new)]), rownames(clinical.values.prev))) == 0)


  indivs.with.cbc = c(rownames(clinical.values.prev), rownames(clinical.values.new))
  orig.cbc = read.csv(ORIG.CBC.PATH, stringsAsFactors=F)
  orig.cbc = orig.cbc[1:184,]
  rownames(orig.cbc) = paste0('N', orig.cbc[,1])
  orig.cbc = orig.cbc[,3:22]
  colnames(orig.cbc) = c('WBC', 'RBC', 'Hemoglobin', 'Hematocrit', 'MCV', 'MCH', 'MCHC', 'RDW', 'Platelets', 'MPV',
	                 'Neutro%', 'Lympho%', 'Mono%', 'Eosinophils%', 'Basophils%', 'Neutro#', 'Lympho#', 'Mono#', 'Eosinophils#', 'Basophils#')

  corrected.cbc.all.samples = rbind(rbind(clinical.values.prev, clinical.values.new), orig.cbc[setdiff(rownames(orig.cbc), c(indivs.with.cbc, indiv.id.map[indivs.with.cbc])),])
  corrected.cbc.fil = corrected.cbc.all.samples[rowSums(is.na(corrected.cbc.all.samples)) < 20,]
  should.change = rownames(corrected.cbc.fil) %in% names(indiv.id.map)
  rownames(corrected.cbc.fil)[should.change] = indiv.id.map[rownames(corrected.cbc.fil)[should.change]]

  return(corrected.cbc.fil)
}

get.indiv.id.map <- function() {
  indiv.id.map.df = read.table(INDIV.ID.MAP.PATH)
  indiv.id.map = indiv.id.map.df$old_id
  names(indiv.id.map) = indiv.id.map.df$new_id
  return(indiv.id.map)
}

fig2.rdw.vs.mut <- function(fig.dir=BASE.FIG.DIR) {
  cohort = read.csv(file.path(RDW.MUT.DIR, 'cohort.csv'), stringsAsFactors=F)
  controls = cohort$Patient_ID[cohort$RDW == 0]
  rdws = cohort$Patient_ID[cohort$RDW == 1]

  arch.muts.df = read.xlsx(file.path(RDW.MUT.DIR, 'arch.xlsx'), 'ARCH mutations RDW and controls')
  arch.hotspots.df = read.xlsx(file.path(RDW.MUT.DIR, 'arch.xlsx'), 'ARCH hotspots RDW and controls')
  amplicon.df = read.xlsx(file.path(RDW.MUT.DIR, 'arch.xlsx'), 'Amplicon RDW and controls')

  arch.muts.df$mut_id = sprintf('%s_%s_%s', arch.muts.df$PATIENT_ID, arch.muts.df$CHR, arch.muts.df$POS)
  arch.hotspots.df$mut_id = sprintf('%s_%s_%s', arch.hotspots.df$PATIENT_ID, arch.hotspots.df$CHR, arch.hotspots.df$POS)
  amplicon.df$mut_id = sprintf('%s_%s_%s', amplicon.df$Patient.ID, amplicon.df$Chr, amplicon.df$Position)

  amplicon.df.fil = data.frame(indiv=amplicon.df$Patient.ID, cohort=amplicon.df$Cohort, chr=amplicon.df$Chr, pos=amplicon.df$Position, 
                               hotspot=amplicon.df$Hotspot, src='amplicon', mut_id=amplicon.df$mut_id, vaf=amplicon.df$VAF, gene=amplicon.df$Gene)
  amplicon.df.fil = filter(amplicon.df.fil, cohort %in% c('RDW', 'Control'))
  arch.muts.df.fil = data.frame(indiv=arch.muts.df$PATIENT_ID, cohort=arch.muts.df$COHORT, chr=arch.muts.df$CHR, pos=arch.muts.df$POS, 
                               hotspot=0, src='arch_mut', mut_id=arch.muts.df$mut_id, vaf=arch.muts.df$VAF, gene=arch.muts.df$GENE)
  stopifnot(all(arch.muts.df.fil$cohort %in% c('RDW', 'Control')))
  arch.hotspots.df.fil = data.frame(indiv=arch.hotspots.df$PATIENT_ID, cohort=arch.hotspots.df$COHORT, chr=arch.hotspots.df$CHR, pos=arch.hotspots.df$POS, 
                               hotspot=1, src='arch_hotspot', mut_id=arch.hotspots.df$mut_id, vaf=arch.hotspots.df$VAF, gene=arch.hotspots.df$GENE)
  stopifnot(all(arch.hotspots.df.fil$cohort %in% c('RDW', 'Control')))

  all.muts.df = rbind(amplicon.df.fil, arch.muts.df.fil[!(arch.muts.df.fil$mut_id %in% amplicon.df.fil$mut_id),], 
                                       arch.hotspots.df.fil[!(arch.hotspots.df.fil$mut_id %in% amplicon.df.fil$mut_id),])

  mut.id.to.vaf = tapply(1:nrow(all.muts.df), all.muts.df$mut_id, function(indices) median(all.muts.df[indices, 'vaf']))
  all.muts.df$med_vaf = mut.id.to.vaf[all.muts.df$mut_id]
  all.muts.df.fil = all.muts.df[!duplicated(all.muts.df$mut_id),]

  # fig 2F
  mut.rdws = unique(all.muts.df.fil$indiv[all.muts.df.fil$cohort == 'RDW'])
  mut.controls = unique(all.muts.df.fil$indiv[all.muts.df.fil$cohort == 'Control'])
  indiv.gene.id = paste0(all.muts.df.fil$indiv, '_', all.muts.df.fil$gene)
  num.muts.per.cohort = table(all.muts.df.fil[!duplicated(indiv.gene.id), 'cohort'], all.muts.df.fil[!duplicated(indiv.gene.id), 'gene'])
  num.muts.per.cohort.sorted = num.muts.per.cohort[,order(colSums(num.muts.per.cohort), decreasing=T)]
  png(file.path(fig.dir, 'arch_rdw.png'), height=600, width=700)
  barplot(num.muts.per.cohort.sorted / length(controls), beside=T, las=2, col=1:2, width=1.5)
  dev.off()

  # tests
  rdw.mut.tbl = table(c(controls, rdws) %in% rdws, c(controls, rdws) %in% c(mut.rdws, mut.controls))
  print(fisher.test(rdw.mut.tbl))

}

analyze.n122 <- function(our.assign, our.annotation, our.legc, fig.dir=BASE.FIG.DIR) {
  cur.sample.num = 'N122'
  metadata.file.path = N122.GOT.PATH
  load(metadata.file.path, v=T)
  cur.exp.name = 'demux_11_01_21_2'
  cur.exp.cells = rownames(our.cdata$obs)[our.cdata$obs$exp_name == cur.exp.name]
  cur.exp.indiv.cells = cur.exp.cells[our.cdata$obs[cur.exp.cells, 'indiv_id'] == cur.sample.num]
  rownames(metadata) = sprintf('%s_%s', cur.exp.name, rownames(metadata))
  print(table(as.character(our.cdata$obs[cur.exp.cells, 'indiv_id']), metadata[cur.exp.cells, 'Genotype_Approx_Match_No_Gene_Thresh']))

  mut.indiv.cells = intersect(rownames(metadata)[metadata$Genotype_Approx_Match_No_Gene_Thresh == 'MUT'], cur.exp.indiv.cells)
  cur.used.ctypes.ord = c('purple', 'red', 'plum', '#EEBB6E', "#7F9D00", 'gold', 'brown', "#6E1C71", 'blue', 'darkblue', 'lightblue')
  mcs.for.comparison = names(our.annotation)[our.annotation %in% cur.used.ctypes.ord]

  wt.indiv.cells = setdiff(intersect(names(our.assign)[our.assign %in% mcs.for.comparison], cur.exp.indiv.cells), mut.indiv.cells)
  #wt.indiv.cells = intersect(rownames(metadata)[metadata$Genotype_Approx_Match_No_Gene_Thresh == 'WT'], cur.exp.indiv.cells)
  mut.indiv.cells.fil = intersect(mut.indiv.cells, names(our.assign)[our.assign %in% mcs.for.comparison])
  wt.indiv.cells.fil = intersect(wt.indiv.cells, names(our.assign)[our.assign %in% mcs.for.comparison])
  our.annotation.factor = as.factor(our.annotation)
  #tmp.cell.types.ord = c('purple', 'red', 'plum', '#EEBB6E', '#7F9D00', 'gold', 'brown', '#6E1C71', 'blue', 'darkblue', '#566CF2', 'lightblue')
  wt.mut.fil.tbl = rbind(table(our.annotation.factor[our.assign[mut.indiv.cells.fil]]), table(our.annotation.factor[our.assign[wt.indiv.cells.fil]]))[,cur.used.ctypes.ord]
  wt.mut.fil.tbl.norm = wt.mut.fil.tbl / rowSums(wt.mut.fil.tbl)

  # edf 8B
  tmp.clp.types = c('blue', 'darkblue', 'lightblue', '#6E1C71')
  tmp.non.clp.types = setdiff(cur.used.ctypes.ord, tmp.clp.types)
  print(fisher.test(rbind(rowSums(wt.mut.fil.tbl[,tmp.clp.types]), rowSums(wt.mut.fil.tbl[,tmp.non.clp.types]))))
  png(file.path(fig.dir, 'barplot_122.png'))
  barplot(t(wt.mut.fil.tbl.norm[2:1,]), col=colnames(wt.mut.fil.tbl.norm))
  dev.off()
}

sync.analysis <- function(our.annotation, our.assign, our.cdata, our.legc, our.cdata.x, indiv.to.used.samp, cell.to.meas, cell.to.samp, samp.ctype.tbl) {
  strict.traj.mcs = names(our.annotation)[our.annotation %in% c('gold', 'plum', '#EEBB6E', 'brown')]
  traj.cells = names(our.assign)[our.assign %in% strict.traj.mcs]
  avp.cors = apply(our.legc[,strict.traj.mcs], 1, cor, our.legc['AVP', strict.traj.mcs])
  gata1.cors = apply(our.legc[,strict.traj.mcs], 1, cor, our.legc['GATA1', strict.traj.mcs])
  avp.genes.fil = setdiff(names(avp.cors)[avp.cors > 0.6], NA)
  gata1.genes.fil = setdiff(names(gata1.cors)[gata1.cors > 0.7], NA)

  avp.genes.fil = avp.genes.fil[!our.cdata$var[avp.genes.fil, 'lateral_gene']]
  gata1.genes.fil = gata1.genes.fil[!our.cdata$var[gata1.genes.fil, 'lateral_gene']]
  avp.genes.fil = avp.genes.fil[rowMeans(our.legc[avp.genes.fil,strict.traj.mcs]) < -10]
  gata1.genes.fil = gata1.genes.fil[rowMeans(our.legc[gata1.genes.fil,strict.traj.mcs]) < -10]
  stopifnot(length(intersect(avp.genes.fil, gata1.genes.fil)) == 0)

  sync.module.df = data.frame(gene=names(avp.cors), avp_cor=avp.cors, gata1_cor=gata1.cors, 
                              is_lateral=our.cdata$var[names(avp.cors), 'lateral_gene'], 
			      mebemp_exp=rowMeans(our.legc[names(avp.cors), strict.traj.mcs]), 
			      is_avp_sig=names(avp.cors) %in% avp.genes.fil, 
			      is_gata1_sig=names(gata1.cors) %in% gata1.genes.fil)
  write.csv(sync.module.df, file=file.path(SUPP.TABLE.DIR, 's8_sync_genes.csv'), quote=F, row.names=F)

  avp.cell.cov = rowSums(our.cdata.x[traj.cells, avp.genes.fil])
  gata1.cell.cov = rowSums(our.cdata.x[traj.cells, gata1.genes.fil])
  total.cell.cov = rowSums(our.cdata.x[traj.cells,])
  avp.cell.scores = avp.cell.cov / total.cell.cov
  gata1.cell.scores = gata1.cell.cov / total.cell.cov

  avp.bin = sort(rep(1:5, length.out=length(avp.cell.scores)))
  names(avp.bin) = names(avp.cell.scores)[order(avp.cell.scores)]
  avp.bin = avp.bin[traj.cells]
  gata1.bin = sort(rep(1:5, length.out=length(gata1.cell.scores)))
  names(gata1.bin) = names(gata1.cell.scores)[order(gata1.cell.scores)]
  gata1.bin = gata1.bin[traj.cells]
  sync.scores = sapply(indiv.to.used.samp, function(cur.samp) {
    cur.cells = intersect(traj.cells, names(cell.to.samp)[cell.to.samp == cur.samp])
    high.gata1.cells = cur.cells[gata1.bin[cur.cells] > 3]
    mean(avp.bin[high.gata1.cells] > 2)
  })
  names(sync.scores) = indiv.to.used.samp

  gata1.bin.hires = sort(rep(1:20, length.out=length(gata1.cell.scores)))
  names(gata1.bin.hires) = names(gata1.cell.scores)[order(gata1.cell.scores)]
  gata1.bin.hires = as.factor(gata1.bin.hires[traj.cells])
  #avp.bin.hires = sort(rep(1:10, length.out=length(avp.cell.scores)))
  avp.bin.hires = sort(rep(1:20, length.out=length(avp.cell.scores)))
  names(avp.bin.hires) = names(avp.cell.scores)[order(avp.cell.scores)]
  avp.bin.hires = as.factor(avp.bin.hires[traj.cells])

  cells.to.plot = traj.cells[(avp.cell.cov[traj.cells] + gata1.cell.cov[traj.cells]) > 10]
  indivs.for.sync.plot = c('N122', 'N172', 'N16', 'N86', 'N98', 'N121')

  smooth.matrix <- function(mat) {
    tmp.mat = matrix(NA, ncol=ncol(mat) + 2, nrow=nrow(mat) + 2)
    tmp.mat[2:(nrow(mat) + 1), 2:(ncol(mat) + 1)] = mat
    ret.mat = matrix(NA, ncol=ncol(mat), nrow=nrow(mat))
    for (i in 1:nrow(ret.mat)) {
      for (j in 1:ncol(ret.mat)) {
        ret.mat[i, j] = mean(tmp.mat[i:(i+2), j:(j+2)], na.rm=T)
      }
    }
    return(ret.mat)
  }

  # fig 3I
  set.seed(42)
  sync.fig.dir = file.path(fig.dir, 'sync_indiv_figs')
  dir.create(sync.fig.dir, showWarnings=F)
  for (cur.indiv in indivs.for.sync.plot) {
    print(cur.indiv)
    cells.to.plot2 = cells.to.plot[cell.to.samp[cells.to.plot] == indiv.to.used.samp[cur.indiv]]
    cur.filename = file.path(sync.fig.dir, paste0(cur.indiv, '.png'))
    tmp.tbl = 0.5 + table(avp.bin.hires[cells.to.plot2], gata1.bin.hires[cells.to.plot2])[20:1,]
    tmp.tbl.norm = tmp.tbl / sum(tmp.tbl)
    pheatmap(pmax(pmin(log2(smooth.matrix(tmp.tbl.norm)), -7), -11), cluster_rows=F, cluster_cols=F, filename=cur.filename, height=5, width=6, breaks=seq(-11, -7, length.out=100))
  }

  used.samps = intersect(rownames(samp.ctype.tbl)[rowSums(samp.ctype.tbl) >= 500], indiv.to.used.samp)
  samp.to.indiv = sapply(strsplit(all.samp, '-'), function(x) x[2])
  names(samp.to.indiv) = all.samp

  used.males = samp.to.indiv[used.samps][is.male[samp.to.indiv[used.samps]]]
  used.male.samps = used.samps[is.male[samp.to.indiv[used.samps]]]
  used.males = samp.to.indiv[used.male.samps]

  indiv.id.map = get.indiv.id.map()
  cbc.new = get.clinical.values.new(indiv.id.map)

  print(kruskal.test(cbc.new[used.males,"MCV"], g=sync.scores[used.male.samps]> .2368))

  # fig 3J
  png(file.path(fig.dir, 'sync_rbc_mcv.png'), width=500, height=500)
  plot(cbc.new[used.males,"RBC"]*10, cbc.new[used.males,"MCV"], pch=19, col=ifelse(sync.scores[used.male.samps]> .2368,"red","black"),cex=2)
  dev.off()

  # fig 3K
  samp.ctype.tbl.norm = samp.ctype.tbl / rowSums(samp.ctype.tbl)
  png(file.path(fig.dir, 'sync_vs_cell_type_freq.png'), height=750, width=1000)
  temp.ret = permut_test(used.male.samps, sync.scores[used.male.samps], 'sync', samp.ctype.tbl.norm, win=3, num.perm=1e4)
  dev.off()

  sync.df = data.frame(indiv_id=samp.to.indiv[used.samps], sync_score=sync.scores[used.samps])
  write.csv(sync.df, file=file.path(SUPP.TABLE.DIR, 'sync_scores.csv'), quote=F, row.names=F)

  # replicates
  sync.scores.samps = sapply(c(unlist(bio.reps[,1]), unlist(bio.reps[,2])), function(cur.samp) {
    cur.cells = intersect(traj.cells, names(cell.to.samp)[cell.to.samp == cur.samp])
    high.gata1.cells = cur.cells[gata1.bin[cur.cells] > 3]
    mean(avp.bin[high.gata1.cells] > 2)
  })
  sync.scores.meas = sapply(c(unlist(tech.reps[,1]), unlist(tech.reps[,2])), function(cur.meas) {
    cur.cells = intersect(traj.cells, names(cell.to.meas)[cell.to.meas == cur.meas])
    high.gata1.cells = cur.cells[gata1.bin[cur.cells] > 3]
    mean(avp.bin[high.gata1.cells] > 2)
  })

  cur.lim = range(c(sync.scores.samps[unlist(bio.reps[,1])], sync.scores.samps[unlist(bio.reps[,2])], 
                    sync.scores.meas[unlist(tech.reps[,1])], sync.scores.meas[unlist(tech.reps[,2])]))
  # edf 8H
  png(file.path(fig.dir, 'sync_bio_reps.png'), width=600, height=600)
  plot(sync.scores.samps[unlist(bio.reps[,1])], sync.scores.samps[unlist(bio.reps[,2])], cex=3, pch=19, xlim=cur.lim, ylim=cur.lim)
  abline(b=1, a=0, col=2, lwd=3)
  dev.off()
  png(file.path(fig.dir, 'sync_tech_reps.png'), width=600, height=600)
  plot(sync.scores.meas[unlist(tech.reps[,1])], sync.scores.meas[unlist(tech.reps[,2])], cex=3, pch=19, xlim=cur.lim, ylim=cur.lim)
  abline(b=1, a=0, col=2, lwd=3)
  dev.off()

}

get.meas.de.expression <- function(our.cdata, our.cdata.x, our.annotation, our.assign, cell.to.meas, tmp.measurements, binned.pooled, bin.per.mc, cur.ctypes, cur.thresh, new.batch.coeffs) {

  cur.traj.mcs = names(our.annotation)[our.annotation %in% cur.ctypes]
  cur.cells = names(our.assign)[our.assign %in% cur.traj.mcs]

  # new code for pooling
  #tmp.measurements = unique(cell.to.meas[cur.cells])
  all.sel.cells = list()
  cur.exp.pooled = NULL
  set.seed(42)
  for (tmp.meas in tmp.measurements) {
    print(tmp.meas)
    print(Sys.time())
    cells = intersect(names(cell.to.meas)[cell.to.meas == tmp.meas], cur.cells)
    cur.meas.mat = our.cdata.x[cells,,drop=F]
    
    cur.num.ds = max(quantile(rowSums(cur.meas.mat), 0.05), 500)
    #cur.num.ds = 500
    if (sum(rowSums(cur.meas.mat) > cur.num.ds) < 3) next
    tmp.ds = scm_downsamp(t(cur.meas.mat), cur.num.ds) 
    cells.fil = colnames(tmp.ds)
    #print(length(cells.fil))
    all.sel.cells[[tmp.meas]] = cells.fil
    cur.exp.pooled = cbind(cur.exp.pooled, rowSums(tmp.ds[,cells.fil]))
  }
  colnames(cur.exp.pooled) = names(all.sel.cells)
  cur.exp.pooled.norm = t(t(cur.exp.pooled) / colSums(cur.exp.pooled))
  cur.meass.legc.raw = log2(1e-5 + cur.exp.pooled.norm)
  cur.meass = names(which(colSums(cur.exp.pooled) > cur.thresh))
  cur.meass.legc = cur.meass.legc.raw

  cur.exp.pooled.traj.norm = do.call(cbind, lapply(cur.meass, function(cur.meas) {
    cur.meas.cells = all.sel.cells[[cur.meas]]
    tmp.pool = rowMeans(binned.pooled[,bin.per.mc[as.character(our.assign[cur.meas.cells])]])
    expected.legc = log2(1e-5 + tmp.pool / sum(tmp.pool))
    observed.legc = cur.meass.legc[,cur.meas]
    return(observed.legc - expected.legc)
  }))
  colnames(cur.exp.pooled.traj.norm) = cur.meass

  meas.src.tbl = table(cell.to.meas, as.character(our.cdata$obs$src))
  stopifnot(all(rowSums(meas.src.tbl > 0) == 1))
  meas.to.src = colnames(meas.src.tbl)[apply(meas.src.tbl, 1, which.max)]
  names(meas.to.src) = rownames(meas.src.tbl)

  meas.to.is.late.exp = meas.to.src[cur.meass] == 'new'
  meas.to.is.late.exp[startsWith(names(meas.to.is.late.exp), 'demux_22_02_21_1')] = T
  meas.to.is.late.exp[startsWith(names(meas.to.is.late.exp), 'demux_22_02_21_2')] = T
  cur.late.meass = names(meas.to.is.late.exp)[meas.to.is.late.exp]
  cur.early.meass = names(meas.to.is.late.exp)[!meas.to.is.late.exp]
  cur.exp.pooled.traj.norm.fixed = cur.exp.pooled.traj.norm[, cur.meass]
  cur.exp.pooled.traj.norm.fixed[,cur.late.meass] = cur.exp.pooled.traj.norm.fixed[,cur.late.meass] - new.batch.coeffs
  return(cur.exp.pooled.traj.norm.fixed)
}

de.analysis <- function(our.cdata, our.cdata.x, our.mdata, our.legc, our.annotation, our.assign, cell.to.samp, indiv.to.used.samp, indiv.to.used.samp.fil, fig.dir=BASE.FIG.DIR) {
  all.exp.pooled = list()
  all.exp.pooled.traj.norm = list()
  all.exp.pooled.traj.norm.fixed = list()
  all.exp.pooled.traj.norm.unfil.fixed = list()
  all.meas.traj.norm.fixed = list()
  all.samps.legcs = list()
  all.is.bio.rep.used = list()

  samp.src.tbl = table(cell.to.samp, as.character(our.cdata$obs$src))
  stopifnot(all(rowSums(samp.src.tbl > 0) == 1))
  samp.to.src = colnames(samp.src.tbl)[apply(samp.src.tbl, 1, which.max)]
  names(samp.to.src) = rownames(samp.src.tbl)

  tech.reps = get.tech.replicates(cell.to.meas, names(cell.to.meas))

  for (i in 1:2) {
    #cur.ctypes = list(c('gold', 'plum', '#EEBB6E', 'brown'), c('blue', 'darkblue', '#6E1C71'))[[i]]
    cur.ctypes = list('gold', 'blue')[[i]]
    cur.traj.mcs = names(our.annotation)[our.annotation %in% cur.ctypes]
    cur.cells = names(our.assign)[our.assign %in% cur.traj.mcs]

    # new code for pooling
    tmp.samps = unique(cell.to.samp[cur.cells])
    all.sel.cells = list()
    cur.exp.pooled = NULL
    set.seed(42)
    for (tmp.samp in tmp.samps) {
      print(tmp.samp)
      print(Sys.time())
      cells = intersect(names(cell.to.samp)[cell.to.samp == tmp.samp], cur.cells)
      cur.samp.mat = our.cdata.x[cells,,drop=F]
      
      cur.num.ds = max(quantile(rowSums(cur.samp.mat), 0.05), 500)
      #cur.num.ds = 500
      #cur.num.ds = max(quantile(rowSums(cur.samp.mat), 0.05), 2000)
      if (sum(rowSums(cur.samp.mat) > cur.num.ds) < 3) next
      tmp.ds = scm_downsamp(t(cur.samp.mat), cur.num.ds) 
      cells.fil = colnames(tmp.ds)
      #print(length(cells.fil))
      all.sel.cells[[tmp.samp]] = cells.fil
      cur.exp.pooled = cbind(cur.exp.pooled, rowSums(tmp.ds[,cells.fil]))
    }
    colnames(cur.exp.pooled) = names(all.sel.cells)
    all.exp.pooled[[i]] = cur.exp.pooled
    #cur.thresh = c(3e5, 1e5)[i]
    cur.thresh = c(2e5, 1e5)[i]
    cur.exp.pooled.norm = t(t(cur.exp.pooled) / colSums(cur.exp.pooled))
    #cur.samps.legc = log2(1e-5 + cur.exp.pooled.norm)
    cur.samps.legc.raw = log2(1e-5 + cur.exp.pooled.norm)
    cur.samps.legc = cur.samps.legc.raw
    #ds.thresh=0.2
    #rs.thresh=0.28
    cur.samps = intersect(names(which(colSums(cur.exp.pooled) > cur.thresh)), indiv.to.used.samp)
    cur.indivs = samp.to.indiv[cur.samps]

    cur.samps.legc.fil = cur.samps.legc[,cur.samps]
    all.samps.legcs[[i]] = cur.samps.legc.fil


    if (i == 1) {
      mebemp.avp.scores.hi.res = our.legc['AVP', our.assign[cur.cells]]
      names(mebemp.avp.scores.hi.res) = cur.cells
      num.hi.res.bins = 30
      mebemp.bin.hi.res = sort(rep(1:num.hi.res.bins, length.out=length(mebemp.avp.scores.hi.res)))
      names(mebemp.bin.hi.res) = names(mebemp.avp.scores.hi.res)[order(mebemp.avp.scores.hi.res, decreasing=T)]
      cur.bin.hi.res = mebemp.bin.hi.res
    } else {
      clp.dntt.scores.hi.res = our.legc['DNTT', our.assign[cur.cells]] + our.legc['VPREB1', our.assign[cur.cells]]
      names(clp.dntt.scores.hi.res) = cur.cells
      num.hi.res.bins = 6
      clp.bin.hi.res = sort(rep(1:num.hi.res.bins, length.out=length(clp.dntt.scores.hi.res)))
      names(clp.bin.hi.res) = names(clp.dntt.scores.hi.res)[order(clp.dntt.scores.hi.res, decreasing=T)]
      cur.bin.hi.res = clp.bin.hi.res
    }
    tmp.tbl = table(our.assign[cur.cells], cur.bin.hi.res[cur.cells])
    bin.per.mc = colnames(tmp.tbl)[apply(tmp.tbl, 1, which.max)]
    names(bin.per.mc) = rownames(tmp.tbl)

    set.seed(42)
    our.mdata.x.ds = t(scm_downsamp(t(our.mdata.x), 9e4))
    stopifnot(nrow(our.mdata.x.ds) == nrow(our.mdata.x))
    binned.pooled = do.call(cbind, tapply(names(bin.per.mc), bin.per.mc, function(cur.mcs) {
      colMeans(our.mdata.x.ds[cur.mcs,])
    }))
    binned.pooled = binned.pooled[,as.character(1:num.hi.res.bins)]
 
    # new traj norm code
    cur.exp.pooled.traj.norm.unfil = do.call(cbind, lapply(colnames(cur.samps.legc), function(cur.samp) {
      cur.samp.cells = all.sel.cells[[cur.samp]]
      #expected.legc = rowMeans(binned.legc[,bin.per.mc[as.character(our.assign[cur.samp.cells])]])
      tmp.pool = rowMeans(binned.pooled[,bin.per.mc[as.character(our.assign[cur.samp.cells])]])
      expected.legc = log2(1e-5 + tmp.pool / sum(tmp.pool))
      observed.legc = cur.samps.legc[,cur.samp]
      return(observed.legc - expected.legc)
    }))
    colnames(cur.exp.pooled.traj.norm.unfil) = colnames(cur.samps.legc)
    cur.exp.pooled.traj.norm = cur.exp.pooled.traj.norm.unfil[,cur.samps]
    all.exp.pooled.traj.norm[[i]] = cur.exp.pooled.traj.norm

    cur.kruskal.pvals.df = read.table(BATCH.EFFECT.PVALS.PATH, header=T)
    cur.kruskal.pvals = cur.kruskal.pvals.df[,2]
    names(cur.kruskal.pvals) = cur.kruskal.pvals.df[,1]

    # also filter genes that are too highly correlated to the traj
    if (i == 1) {
      avp.cors = apply(our.legc[,cur.traj.mcs], 1, cor, -our.legc['AVP', cur.traj.mcs])
      traj.cor.genes = sort(names(avp.cors)[abs(avp.cors) > 0.65])
    } else {
      dntt.cors = apply(our.legc[,cur.traj.mcs], 1, cor, (our.legc['DNTT', cur.traj.mcs] + our.legc['VPREB1', cur.traj.mcs]))
      traj.cor.genes = sort(names(dntt.cors)[abs(dntt.cors) > 0.65])
    }

    # calculate correction factor for new samples
    samp.to.is.late.exp = samp.to.src[cur.samps] == 'new'
    samp.to.is.late.exp[startsWith(names(samp.to.is.late.exp), 'demux_22_02_21_ultima')] = T
    cur.late.samps = names(samp.to.is.late.exp)[samp.to.is.late.exp]
    cur.early.samps = names(samp.to.is.late.exp)[!samp.to.is.late.exp]
    new.batch.coeffs = sapply(1:nrow(cur.exp.pooled.traj.norm), function(j) {
      if (j %% 100 == 0) print(j)
      cur.gene = rownames(cur.exp.pooled.traj.norm)[j]
      cur.df = data.frame(vals=cur.exp.pooled.traj.norm[cur.gene, cur.samps], is_male=sexes[samp.to.indiv[cur.samps]] == 'male', age=ages[samp.to.indiv[cur.samps]], is_new=samp.to.is.late.exp[cur.samps])
      cur.lm.ret = lm(vals ~ is_male + age + is_new, data=cur.df)
      return(cur.lm.ret$coeff[4])
    })
    names(new.batch.coeffs) = rownames(cur.exp.pooled.traj.norm)

    cur.exp.pooled.traj.norm.fixed = cur.exp.pooled.traj.norm[, cur.samps]
    cur.exp.pooled.traj.norm.fixed[,cur.late.samps] = cur.exp.pooled.traj.norm.fixed[,cur.late.samps] - new.batch.coeffs
    all.exp.pooled.traj.norm.fixed[[i]] = cur.exp.pooled.traj.norm.fixed

    unfil.samp.to.is.late.exp = samp.to.src[colnames(cur.exp.pooled.traj.norm.unfil)] == 'new'
    unfil.samp.to.is.late.exp[startsWith(names(unfil.samp.to.is.late.exp), 'demux_22_02_21_ultima')] = T
    unfil.late.samps = names(unfil.samp.to.is.late.exp)[unfil.samp.to.is.late.exp]
    unfil.early.samps = names(unfil.samp.to.is.late.exp)[!unfil.samp.to.is.late.exp]
    cur.exp.pooled.traj.norm.unfil.fixed = cur.exp.pooled.traj.norm.unfil
    cur.exp.pooled.traj.norm.unfil.fixed[,unfil.late.samps] = cur.exp.pooled.traj.norm.unfil.fixed[,unfil.late.samps] - new.batch.coeffs

    all.exp.pooled.traj.norm.unfil.fixed[[i]] = cur.exp.pooled.traj.norm.unfil.fixed                                                                                                                    


    late.mod.genes = c('PEX6', 'NARFL', 'CNN2', 'CTSD', 'LAT2', 'RNF187', 'ASMTL', 'CBX6', 'ACAP1', 'LTBP4', 'PPP1R15A', 'ABHD17A', 'C19orf48', 'ESRRA', 'KCNN4', 'PQLC1', 'SLC25A29', 'ARL16', 'RRP7A', 'CAPG', 'FOSB', 'KLF6', 'RRBP1', 'ZNF83', 'ZNF141', 'PBX2', 'MFGE8', 'GIGYF1', 'ERF', 'ATF6B', 'CBFA2T3', 'FLNA', 'BTBD2', 'PIM1', 'SIL1', 'EPN1')
    late.mod.cors = apply(cur.exp.pooled.traj.norm, 1, cor, colMeans(cur.exp.pooled.traj.norm[late.mod.genes,]))
    late.mod.cor.genes = names(which(abs(late.mod.cors) > 0.5))
    genes.for.analysis = setdiff(rownames(our.legc), c(sort(names(cur.kruskal.pvals)[cur.kruskal.pvals < 0.001]), c(traj.cor.genes, late.mod.cor.genes)))
    #genes.for.analysis = setdiff(rownames(our.legc), c(sort(names(cur.kruskal.pvals)[cur.kruskal.pvals < 0.001]), late.mod.cor.genes))

    # clock 
    tmp.max.exp = apply(cur.samps.legc.fil, 1, max)
    min.exp.for.clock = c(-14.5, -15.5)[i]
    genes.for.clock = intersect(names(tmp.max.exp)[tmp.max.exp > min.exp.for.clock], genes.for.analysis)
    cur.exp.pooled.traj.norm.fixed2 = cur.exp.pooled.traj.norm.fixed[genes.for.clock,]
    cv_model <- cv.glmnet(t(cur.exp.pooled.traj.norm.fixed2), ages[cur.indivs], alpha = 1, nfolds=length(cur.indivs))
    best_lambda <- cv_model$lambda.min
    best_model <- glmnet(t(cur.exp.pooled.traj.norm.fixed2), ages[cur.indivs], alpha = 1, nfolds=length(cur.indivs), lambda=best_lambda)
    tmp.preds = sapply(1:length(cur.indivs), function(j) {
      cur_best_model <- glmnet(t(cur.exp.pooled.traj.norm.fixed2[,-j]), ages[cur.indivs][-j], alpha = 1, lambda=best_lambda)
      best_model_pred = predict(cur_best_model, s = best_lambda, newx = t(cur.exp.pooled.traj.norm.fixed2[,j]))[,1]
      #best_model_pred = predict(best_model, s = best_lambda, newx = t(cur.exp.pooled.traj.norm.fixed2[,j]))[,1]
      return(best_model_pred)
    })
    png(file.path(fig.dir, c('age_norm_exp_clock_unpooled_loocv_mebemp.png', 'age_norm_exp_clock_unpooled_loocv_clp.png')[i]))
    tmp.lim = range(c(ages[cur.indivs], tmp.preds))
    #tmp.cor = cor(ages[cur.indivs], tmp.preds, method='spearman')
    tmp.cor = cor(ages[cur.indivs], tmp.preds)
    plot(ages[cur.indivs], tmp.preds, xlim=tmp.lim, ylim=tmp.lim, pch=19, main=paste('R^2 is:', round(tmp.cor, 3) ** 2), cex=2)
    abline(b=1, a=0, col=2, lwd=3)
    dev.off()

    # clock nested cross validation
    tmp.preds.nested = sapply(1:length(cur.indivs), function(j) {
      #if (j %% 10 == 1) print(j)
      print(j)
      print(Sys.time())
      cur_cv_model <- cv.glmnet(t(cur.exp.pooled.traj.norm.fixed2[,-j]), ages[cur.indivs][-j], alpha = 1, nfolds=length(cur.indivs) - 1)
      cur_best_lambda = cur_cv_model$lambda.min
      cur_best_model <- glmnet(t(cur.exp.pooled.traj.norm.fixed2[,-j]), ages[cur.indivs][-j], alpha = 1, lambda=cur_best_lambda)
      best_model_pred = predict(cur_best_model, s = cur_best_lambda, newx = t(cur.exp.pooled.traj.norm.fixed2[,j]))[,1]
      #best_model_pred = predict(best_model, s = best_lambda, newx = t(cur.exp.pooled.traj.norm.fixed2[,j]))[,1]
      return(best_model_pred)
    })
    # fig 3E, edf 8C
    png(file.path(fig.dir, c('nested_age_norm_exp_clock_unpooled_loocv_mebemp.png', 'nested_age_norm_exp_clock_unpooled_loocv_clp.png')[i]))
    tmp.lim.nested = range(c(ages[cur.indivs], tmp.preds.nested))
    #tmp.cor = cor(ages[cur.indivs], tmp.preds, method='spearman')
    tmp.cor.nested = cor(ages[cur.indivs], tmp.preds.nested)
    plot(ages[cur.indivs], tmp.preds.nested, xlim=tmp.lim.nested, ylim=tmp.lim.nested, pch=19, main=paste('R^2 is:', round(tmp.cor.nested, 3) ** 2), cex=2)
    abline(b=1, a=0, col=2, lwd=3)
    dev.off()



    # find high var genes
    gene.means = sort(rowMeans(cur.exp.pooled.norm[genes.for.analysis, cur.samps]))
    gene.sds = apply(cur.exp.pooled.traj.norm, 1, sd)[names(gene.means)]
    gene.sds.smoothed = rollmean(gene.sds, k=100, fill='extend')
    #plot(log2(1e-5 + gene.means), gene.sds, pch=19, cex=0.5)
    #grid()
    #lines(log2(1e-5 + gene.means), gene.sds.smoothed + 0.08, col=2)
    #abline(v=-8, col=2)
    high.var.genes = names(gene.sds)[gene.sds > gene.sds.smoothed + 0.08 & log2(1e-5 + gene.means) < -8]
    #print(sort(high.var.genes))
    #cur.gene.cor = tgs_cor(data.matrix(t(cur.indiv.legc.[high.var.genes,])), spearman=T)
    cur.gene.cor = tgs_cor(data.matrix(t(cur.exp.pooled.traj.norm.fixed[high.var.genes,])), spearman=T)
    cur.gene.cor.early = tgs_cor(data.matrix(t(cur.exp.pooled.traj.norm.fixed[high.var.genes, cur.early.samps])), spearman=T)
    shades = colorRampPalette(c("darkblue", "blue", "white","red", "darkred"))(100)
    #shades[25:60] = 'white'
    shades[40:60] = 'white'
    diag(cur.gene.cor) = NA
    pheatmap.ret = pheatmap(cur.gene.cor, breaks=seq(-1, 1, length.out=100), color=shades)
    num.cls = 30
    tmp.cls = cutree(pheatmap.ret$tree_col, num.cls)
    #mean.cls.cor = sapply(1:num.cls, function(cur.cls) mean(cur.gene.cor[tmp.cls == cur.cls, tmp.cls == cur.cls], na.rm=T))
    #high.var.genes2 = setdiff(names(tmp.cls)[mean.cls.cor[tmp.cls] > 0.2], names(tmp.cls)[tmp.cls %in% tmp.cls[c('BIRC3', 'PCDH9')]])

    high.var.genes2 = high.var.genes[sapply(high.var.genes, function(cur.gene) mean(cur.gene.cor[cur.gene, tmp.cls == tmp.cls[cur.gene]], na.rm=T)) > 0.2]
    mean.cls.cor = tapply(high.var.genes2, tmp.cls[high.var.genes2], function(cur.genes) mean(cur.gene.cor[cur.genes, cur.genes], na.rm=T))
    mean.cls.cor.early = tapply(high.var.genes2, tmp.cls[high.var.genes2], function(cur.genes) mean(cur.gene.cor.early[cur.genes, cur.genes], na.rm=T))
    high.var.genes2.fil = intersect(high.var.genes2[mean.cls.cor[as.character(tmp.cls[high.var.genes2])] > 0.25], 
                                high.var.genes2[mean.cls.cor.early[as.character(tmp.cls[high.var.genes2])] > 0.25])

    genes.per.cls = table(tmp.cls[high.var.genes2.fil])
    high.var.genes3 = high.var.genes2.fil[genes.per.cls[as.character(tmp.cls[high.var.genes2.fil])] > 2]
    cur.gene.cor.trim = cur.gene.cor[high.var.genes3, high.var.genes3]

    # all high var genes
    pheatmap.ret = pheatmap(cur.gene.cor, breaks=seq(-1, 1, length.out=100), color=shades, border_color=F, fontsize=7)
    pheatmap(cur.gene.cor[pheatmap.ret$tree_row$ord, pheatmap.ret$tree_col$ord], breaks=seq(-1, 1, length.out=100), color=shades, 
                            show_colnames=T, show_rownames=F, cluster_rows=F, cluster_cols=F,
                            border_color=F, fontsize=2, filename=file.path(fig.dir, c('mebemp_modules_heatmap_unfil.png', 'clp_modules_heatmap_unfil.png')[i]), width=10.1, height=5)
    pheatmap(cur.gene.cor.early[pheatmap.ret$tree_row$ord, pheatmap.ret$tree_col$ord], breaks=seq(-1, 1, length.out=100), color=shades, 
                            show_colnames=T, show_rownames=F, cluster_rows=F, cluster_cols=F,
                            border_color=F, fontsize=2, filename=file.path(fig.dir, c('mebemp_modules_heatmap_unfil_early.png', 'clp_modules_heatmap_unfil_early.png')[i]), width=10.1, height=5)

    # "good" high var genes
    pheatmap.ret = pheatmap(cur.gene.cor.trim, breaks=seq(-1, 1, length.out=100), color=shades, border_color=F, fontsize=7)
    pheatmap(cur.gene.cor.trim[pheatmap.ret$tree_row$ord, pheatmap.ret$tree_col$ord], breaks=seq(-1, 1, length.out=100), color=shades, 
                            show_colnames=T, show_rownames=F, cluster_rows=F, cluster_cols=F,
                            border_color=F, fontsize=2, filename=file.path(fig.dir, c('mebemp_modules_heatmap_with_names.png', 'clp_modules_heatmap_with_names.png')[i]), width=10.1, height=5)
    pheatmap(cur.gene.cor.trim[pheatmap.ret$tree_row$ord, pheatmap.ret$tree_col$ord], breaks=seq(-1, 1, length.out=100), color=shades, 
                            show_colnames=F, show_rownames=F, cluster_rows=F, cluster_cols=F,
                            border_color=F, fontsize=7, filename=file.path(fig.dir, c('mebemp_modules_heatmap.png', 'clp_modules_heatmap.png')[i]), width=5.1, height=5)
    pheatmap(cur.gene.cor.early[colnames(cur.gene.cor.trim)[pheatmap.ret$tree_row$ord], rownames(cur.gene.cor.trim)[pheatmap.ret$tree_col$ord]], breaks=seq(-1, 1, length.out=100), color=shades, 
                            show_colnames=F, show_rownames=F, cluster_rows=F, cluster_cols=F,
                            border_color=F, fontsize=7, filename=file.path(fig.dir, c('mebemp_modules_heatmap_early.png', 'clp_modules_heatmap_early.png')[i]), width=5.1, height=5)

    if (i == 1) {
      # write which modules are shown and which are not
      removed.cls = c('ITGA2B', 'PCDH9', 'STAT1', 'CNRIP1', 'DUSP6')
      kept.cls = c('DNTT', 'HIST2H3D', 'UTY', 'XIST', 'CP', 'EMP1',  'CEBPA', 'LMNA', 'MPO')
    } else {
      removed.cls = c('ITGA2B', 'ADAM8', 'BCL9L', 'CDK9', 'PCDH9', 'NLRP1', 'GBP1')
      kept.cls = c('UTY', 'LMNA', 'XIST') 
    }
    remaining.cls = names(table(tmp.cls[high.var.genes3]))
    genes.in.kept.cls = c()
    for (j in seq_along(remaining.cls)) {
      cur.cls.genes = high.var.genes3[tmp.cls[high.var.genes3] == remaining.cls[j]]
      cur.representative = intersect(cur.cls.genes, c(removed.cls, kept.cls))
      stopifnot(length(cur.representative) == 1)
      if (cur.representative %in% kept.cls) {
        genes.in.kept.cls = c(genes.in.kept.cls, cur.cls.genes)
      }
    }
  
    cur.gene.cor.trim2 = cur.gene.cor.trim[genes.in.kept.cls, genes.in.kept.cls]
    pheatmap.ret = pheatmap(cur.gene.cor.trim2, breaks=seq(-1, 1, length.out=100), color=shades, border_color=F, fontsize=7)
    pheatmap(cur.gene.cor.trim2[pheatmap.ret$tree_row$ord, pheatmap.ret$tree_col$ord], breaks=seq(-1, 1, length.out=100), color=shades, 
                            show_colnames=T, show_rownames=F, cluster_rows=F, cluster_cols=F,
                            border_color=F, fontsize=2, filename=file.path(fig.dir, c('mebemp_modules_heatmap_with_names_fil.png', 'clp_modules_heatmap_with_names_fil.png')[i]), width=10.1, height=5)
    # fig 3F
    pheatmap(cur.gene.cor.trim2[pheatmap.ret$tree_row$ord, pheatmap.ret$tree_col$ord], breaks=seq(-1, 1, length.out=100), color=shades, 
                            show_colnames=F, show_rownames=F, cluster_rows=F, cluster_cols=F,
                            border_color=F, fontsize=7, filename=file.path(fig.dir, c('mebemp_modules_heatmap_fil.png', 'clp_modules_heatmap_fil.png')[i]), width=5.1, height=5)
  

    # tech and bio reps
    bio.reps = get.biological.replicates(all.samp, indiv.to.used.samp.fil)
    is.bio.rep.used = colSums(cur.exp.pooled)[unlist(bio.reps[,1])] > c(2e5, 1e5)[i] & colSums(cur.exp.pooled)[unlist(bio.reps[,2])] > c(2e5, 1e5)[i]
    all.is.bio.rep.used[[i]] = is.bio.rep.used

    tmp.measurements = c(unlist(tech.reps[,1]), unlist(tech.reps[,2]))
    meas.traj.norm.fixed = get.meas.de.expression(our.cdata, our.cdata.x, our.annotation, our.assign, 
                                     cell.to.meas, tmp.measurements, binned.pooled, bin.per.mc, cur.ctypes, cur.thresh, new.batch.coeffs)
    all.meas.traj.norm.fixed[[i]] =  meas.traj.norm.fixed

  }

  # LMNA module analysis
  cur.exp.pooled.traj.norm.fixed.mebemp = all.exp.pooled.traj.norm.fixed[[1]]
  cur.exp.pooled.traj.norm.fixed.clp = all.exp.pooled.traj.norm.fixed[[2]]
  lmna.clp.cors = apply(cur.exp.pooled.traj.norm.fixed.clp, 1, cor, cur.exp.pooled.traj.norm.fixed.clp['LMNA',], method='spearman')
  lmna.mebemp.cors = apply(cur.exp.pooled.traj.norm.fixed.mebemp, 1, cor, cur.exp.pooled.traj.norm.fixed.mebemp['LMNA',], method='spearman')
  total.lmna.cor = lmna.mebemp.cors + lmna.clp.cors > 0.7
  png(file.path(fig.dir, 'lmna_mebemp_vs_clp_cors.png'))
  plot(lmna.mebemp.cors, lmna.clp.cors, xlab='MEBEMP LMNA correlation', ylab='CLP LMNA correlation')
  points(lmna.mebemp.cors[total.lmna.cor], lmna.clp.cors[total.lmna.cor], col=2)
  dev.off()
  comb.lmna.mod = sort(names(total.lmna.cor)[total.lmna.cor])
  lmna.mod.fil = setdiff(comb.lmna.mod, c(sort(names(cur.kruskal.pvals)[cur.kruskal.pvals < 0.001]), 'VIM', 'KLF6'))
  # compare scores
  common.samps = intersect(colnames(cur.exp.pooled.traj.norm.fixed.mebemp), colnames(cur.exp.pooled.traj.norm.fixed.clp))
  lmna.clp.scores = colMeans(cur.exp.pooled.traj.norm.fixed.clp[lmna.mod.fil, common.samps])
  lmna.mebemp.scores = colMeans(cur.exp.pooled.traj.norm.fixed.mebemp[lmna.mod.fil, common.samps])
  cols.for.lmna.scatter = ifelse(sexes[names(indiv.to.used.samp)] == 'male', 'blue', 'red')
  # fig 3G
  png(file.path(fig.dir, 'lmna_mebemp_vs_clp_scores.png'))
  tmp.lim = range(c(lmna.mebemp.scores, lmna.clp.scores))
  plot(lmna.mebemp.scores, lmna.clp.scores, pch=21, bg=cols.for.lmna.scatter[samp.to.indiv[common.samps]], cex=2.5, xlim=tmp.lim, ylim=tmp.lim)
  dev.off()

  clp.samps = colnames(cur.exp.pooled.traj.norm.fixed.clp)
  clp.tmp.young.males = clp.samps[ages[samp.to.indiv[clp.samps]] <= 50 & sexes[samp.to.indiv[clp.samps]] == 'male']
  clp.tmp.old.males = clp.samps[ages[samp.to.indiv[clp.samps]] > 60 & sexes[samp.to.indiv[clp.samps]] == 'male']
  clp.tmp.young.females = clp.samps[ages[samp.to.indiv[clp.samps]] <= 50 & sexes[samp.to.indiv[clp.samps]] == 'female']
  clp.tmp.old.females = clp.samps[ages[samp.to.indiv[clp.samps]] > 60 & sexes[samp.to.indiv[clp.samps]] == 'female']

  mebemp.samps = colnames(cur.exp.pooled.traj.norm.fixed.mebemp)
  mebemp.tmp.young.males = mebemp.samps[ages[samp.to.indiv[mebemp.samps]] <= 50 & sexes[samp.to.indiv[mebemp.samps]] == 'male']
  mebemp.tmp.old.males = mebemp.samps[ages[samp.to.indiv[mebemp.samps]] > 60 & sexes[samp.to.indiv[mebemp.samps]] == 'male']
  mebemp.tmp.young.females = mebemp.samps[ages[samp.to.indiv[mebemp.samps]] <= 50 & sexes[samp.to.indiv[mebemp.samps]] == 'female']
  mebemp.tmp.old.females = mebemp.samps[ages[samp.to.indiv[mebemp.samps]] > 60 & sexes[samp.to.indiv[mebemp.samps]] == 'female']


  lmna.mebemp.scores.list = list(lmna.mebemp.scores[c(mebemp.tmp.young.males, mebemp.tmp.young.females)], 
                                 lmna.mebemp.scores[c(mebemp.tmp.old.males, mebemp.tmp.old.females)])
  lmna.clp.scores.list = list(lmna.clp.scores[c(clp.tmp.young.males, clp.tmp.young.females)], lmna.clp.scores[c(clp.tmp.old.males, clp.tmp.old.females)])
  print('kruskal clp')
  print(kruskal.test(lmna.clp.scores.list))
  print('kruskal mebemp')
  print(kruskal.test(lmna.mebemp.scores.list))
  # fig 3H
  png(file.path(fig.dir, 'age_sex_lmna_clp_box.png'), height=800, width=400)
  boxplot(lmna.clp.scores.list, lwd=3)
  dev.off()
  png(file.path(fig.dir, 'age_sex_lmna_mebemp_box.png'), height=800, width=400)
  boxplot(lmna.mebemp.scores.list, lwd=3)
  dev.off()

  cur.used.ctypes = c('purple', 'red', 'plum', '#EEBB6E', "#7F9D00", 'gold', 'brown', "#6E1C71", 'blue', 'darkblue', 'lightblue')
  mcs.for.scatters = names(our.annotation)[our.annotation %in% cur.used.ctypes]
  genes.for.scatters = c('ANXA1', 'TAGLN2', 'AHNAK', 'MYADM', 'TSPAN2', 'VIM')
  # edf 8D
  png(file.path(fig.dir, 'lmna_gg_scatters.png'), height=1500, width=1100)
  par(mfrow=c(3, 2), mar=c(6, 6, 2, 2))
  for (cur.gene in genes.for.scatters) {
    plot(our.legc['LMNA', mcs.for.scatters], our.legc[cur.gene, mcs.for.scatters], ylab=cur.gene, xlab='LMNA', pch=21, bg=our.annotation[mcs.for.scatters], cex=2.5)
    box(lwd=3)
  }
  dev.off()

  # edf 8G
  # bio replicates
  is.bio.rep.used.mebemp = all.is.bio.rep.used[[1]]
  is.bio.rep.used.clp = all.is.bio.rep.used[[2]]
  cur.exp.pooled.traj.norm.unfil.fixed.mebemp = all.exp.pooled.traj.norm.unfil.fixed[[1]]
  cur.exp.pooled.traj.norm.unfil.fixed.clp = all.exp.pooled.traj.norm.unfil.fixed[[2]]
  # mebemp
  mebemp.scores1 = colMeans(cur.exp.pooled.traj.norm.unfil.fixed.mebemp[lmna.mod.fil, unlist(bio.reps[is.bio.rep.used.mebemp, 1])])
  mebemp.scores2 = colMeans(cur.exp.pooled.traj.norm.unfil.fixed.mebemp[lmna.mod.fil, unlist(bio.reps[is.bio.rep.used.mebemp, 2])])
  tmp.lim = range(c(mebemp.scores1, mebemp.scores2))
  png(file.path(fig.dir, 'lmna_mebemp_bio_rep.png'))
  plot(mebemp.scores1, mebemp.scores2, xlim=tmp.lim, ylim=tmp.lim, pch=19, cex=2)
  abline(b=1, a=0, col=2)
  dev.off()
  # clp
  clp.scores1 = colMeans(cur.exp.pooled.traj.norm.unfil.fixed.clp[lmna.mod.fil, unlist(bio.reps[is.bio.rep.used.clp, 1])])
  clp.scores2 = colMeans(cur.exp.pooled.traj.norm.unfil.fixed.clp[lmna.mod.fil, unlist(bio.reps[is.bio.rep.used.clp, 2])])
  tmp.lim = range(c(clp.scores1, clp.scores2))
  png(file.path(fig.dir, 'lmna_clp_bio_rep.png'))
  plot(clp.scores1, clp.scores2, xlim=tmp.lim, ylim=tmp.lim, pch=19, cex=2)
  abline(b=1, a=0, col=2)
  dev.off()

  # technical replicates
  meas.traj.norm.fixed.mebemp = all.meas.traj.norm.fixed[[1]]
  meas.traj.norm.fixed.clp = all.meas.traj.norm.fixed[[2]]
  is.tech.used.mebemp = unlist(tech.reps[,1]) %in% colnames(meas.traj.norm.fixed.mebemp) & unlist(tech.reps[,2]) %in% colnames(meas.traj.norm.fixed.mebemp)
  is.tech.used.clp = unlist(tech.reps[,1]) %in% colnames(meas.traj.norm.fixed.clp) & unlist(tech.reps[,2]) %in% colnames(meas.traj.norm.fixed.clp)
  mebemp.scores1 = colMeans(meas.traj.norm.fixed.mebemp[lmna.mod.fil, unlist(tech.reps[is.tech.used.mebemp, 1])])
  mebemp.scores2 = colMeans(meas.traj.norm.fixed.mebemp[lmna.mod.fil, unlist(tech.reps[is.tech.used.mebemp, 2])])
  clp.scores1 = colMeans(meas.traj.norm.fixed.clp[lmna.mod.fil, unlist(tech.reps[is.tech.used.clp, 1])])
  clp.scores2 = colMeans(meas.traj.norm.fixed.clp[lmna.mod.fil, unlist(tech.reps[is.tech.used.clp, 2])])
  tmp.lim = range(c(mebemp.scores1, mebemp.scores2))
  png(file.path(fig.dir, 'lmna_mebemp_tech_rep.png'))
  plot(mebemp.scores1, mebemp.scores2, xlim=tmp.lim, ylim=tmp.lim, pch=19, cex=2)
  abline(b=1, a=0, col=2)
  dev.off()
  tmp.lim = range(c(clp.scores1, clp.scores2))
  png(file.path(fig.dir, 'lmna_clp_tech_rep.png'))
  plot(clp.scores1, clp.scores2, xlim=tmp.lim, ylim=tmp.lim, pch=19, cex=2)
  abline(b=1, a=0, col=2)
  dev.off()

  # lmna expression across metacells
  mc.lmna.scores = colMeans(our.legc[lmna.mod.fil,])
  #mebemp.mcs = names(our.annotation.amos)[our.annotation.amos %in% c('#EEBB6E', 'gold', 'plum', "#7F9D00", 'red', 'purple')]
  lmna.mebemp.mcs = names(our.annotation)[our.annotation %in% c('#EEBB6E', 'gold', 'plum', "#7F9D00")]
  lmna.brown.mcs = names(our.annotation)[our.annotation == 'brown']
  lmna.clp.mcs = names(our.annotation)[our.annotation %in% c('blue', 'darkblue', '#6E1C71')]
  #tmp.xvals = c(our.legc['AVP', cur.mebemp.mcs], our.legc['AVP', cur.brown.mcs], -(our.legc['AVP', cur.clp.mcs] - max(our.legc['AVP', cur.clp.mcs])) + max(our.legc['AVP', cur.clp.mcs]))
  #plot(tmp.xvals, mc.lmna.scores[names(tmp.xvals)], pch=21, bg=our.annotation[names(tmp.xvals)])

  # edf 8F
  png(file.path(fig.dir, 'avp_vs_lmna_module.png'), width=1100, height=500)
  par(mfrow=c(1, 2))
  cur.ylim = range(mc.lmna.scores[c(lmna.brown.mcs, lmna.mebemp.mcs, lmna.clp.mcs)])
  cur.xlim = range(our.legc['AVP', c(lmna.brown.mcs, lmna.mebemp.mcs, lmna.clp.mcs)])
  #png(file.path(fig.dir, 'avp_vs_lmna_module_mebemp.png'), width=500, height=500)
  plot(our.legc['AVP', c(lmna.brown.mcs, lmna.mebemp.mcs)], mc.lmna.scores[c(lmna.brown.mcs, lmna.mebemp.mcs)], pch=21, bg=our.annotation[c(lmna.brown.mcs, lmna.mebemp.mcs)], xlim=cur.xlim, ylim=cur.ylim, cex=2)
  #dev.off()
  #png(file.path(fig.dir, 'avp_vs_lmna_module_clp.png'), width=500, height=500)
  plot(-our.legc['AVP', c(lmna.brown.mcs, lmna.clp.mcs)], mc.lmna.scores[c(lmna.brown.mcs, lmna.clp.mcs)], pch=21, bg=our.annotation[c(lmna.brown.mcs, lmna.clp.mcs)], xlim=range(-cur.xlim), ylim=cur.ylim, cex=2)
  dev.off()

  mebemp.cells = names(our.assign)[our.assign %in% names(our.annotation)[our.annotation %in% c('gold', '#EEBB6E', 'plum')]]
  mebemp.avp.scores = our.legc['AVP', our.assign[mebemp.cells]]
  names(mebemp.avp.scores) = mebemp.cells
  mebemp.bin = sort(rep(1:10, length.out=length(mebemp.avp.scores)))
  names(mebemp.bin) = names(mebemp.avp.scores)[order(mebemp.avp.scores, decreasing=T)]
  clp.cells = names(our.assign)[our.assign %in% names(our.annotation)[our.annotation %in% c('blue', 'darkblue')]]
  clp.avp.scores = our.legc['AVP', our.assign[clp.cells]]
  names(clp.avp.scores) = clp.cells
  clp.bin = sort(rep(2:4, length.out=length(clp.avp.scores)))
  names(clp.bin) = names(clp.avp.scores)[order(clp.avp.scores, decreasing=T)]
  clp.bin = -clp.bin

  hsc.cells = names(our.assign)[our.assign %in% names(our.annotation)[our.annotation == 'brown']]
  clpe.cells = names(our.assign)[our.assign %in% names(our.annotation)[our.annotation == '#6E1C71']]
  hsc.bin = rep(0, length(hsc.cells))
  names(hsc.bin) = hsc.cells
  clpe.bin = rep(-1, length(clpe.cells))
  names(clpe.bin) = clpe.cells
  main.bin = c(mebemp.bin, hsc.bin, clpe.bin, clp.bin)

  our.cdata.x.lmna = our.cdata.x[,lmna.mod.fil]
  all.traj.lmna.pooled.binned = lapply(1:10, function(i) {
    print(i)
    cur.cells = names(main.bin)[main.bin == as.character(i)]
    cur.cells.fil = cur.cells[cell.to.samp[cur.cells] %in% indiv.to.used.samp]
    cur.exp.pooled = do.call(cbind, tapply(cur.cells.fil, cell.to.samp[cur.cells.fil], function(cells) {colSums(our.cdata.x.lmna[cells,,drop=F])}))
    return(cur.exp.pooled)
  })
  all.cell.cov = rowSums(our.cdata.x)
  total.cov.binned = lapply(1:10, function(i) {
    print(i)
    cur.cells = names(main.bin)[main.bin == as.character(i)]
    cur.cells.fil = cur.cells[cell.to.samp[cur.cells] %in% indiv.to.used.samp]
    tmp = tapply(cur.cells.fil, cell.to.samp[cur.cells.fil], function(cells) {sum(all.cell.cov[cells])})
    cur.total.cov = as.numeric(tmp)
    names(cur.total.cov) = names(tmp)
    return(cur.total.cov)
  })
  samps.for.lmna.heatmap = names(total.cov.binned[[1]])[colSums(do.call(rbind, total.cov.binned)> 1e5) == 10]

  # LMNA score across 10 bins
  lmna.mod.mat = do.call(rbind, lapply(1:10, function(i) rowMeans(log2(1e-5 + t(all.traj.lmna.pooled.binned[[i]]) / total.cov.binned[[i]]))))[,samps.for.lmna.heatmap]
  #lmna.mod.mat = do.call(rbind, all.traj.pooled.lmna.scores)[,sel.indivs]
  lmna.samps.ord = names(sort(colSums(lmna.mod.mat)))
  lmna.indivs.ord = samp.to.indiv[lmna.samps.ord]
  tmp.mat = pmin(pmax(lmna.mod.mat[10:1, lmna.samps.ord], -15.5), -13.5)
  tmp.mat.smoothed = apply(tmp.mat, 2, function(x) c((x[1:(length(x)-1)] + x[2:length(x)]) / 2, x[length(x)]))
  shades = colorRampPalette(c('white', 'red', 'darkred', 'yellow'))(101)
  tmp.df = data.frame(indivs_cols=indivs.cols[lmna.indivs.ord], age=pmin(pmax(ages[lmna.indivs.ord], 20), 85))
  rownames(tmp.df) = lmna.samps.ord
  tmp.cols = unique(indivs.cols)
  names(tmp.cols) = tmp.cols
  age.cols = colorRampPalette(c('white', 'grey', 'black'))(66)
  names(age.cols) = 20:85
  #all.pheatmap = pheatmap(t(tmp.mat.smoothed), breaks=seq(-15.5, -13.5, length.out=100), col=shades, cluster_rows=F, cluster_cols=F, 
  # edf 8E
  all.pheatmap = pheatmap(t(tmp.mat.smoothed), breaks=seq(-15.5, -13.5, length.out=100), col=shades, cluster_rows=F, cluster_cols=F, 
                          annotation_row=tmp.df[colnames(tmp.mat.smoothed),], annotation_colors=list(indivs_cols=tmp.cols, age=age.cols), 
			  filename=file.path(fig.dir, 'lmna_mod_heatmap.png'), height=12, width=10)

  # export LMNA table
  tmp.mebemp.legc = our.legc[,our.annotation[colnames(our.legc)] %in% c('brown', '#EEBB6E', 'gold', 'plum')]
  lmna.mebemp.mc.cors = apply(tmp.mebemp.legc, 1, cor, tmp.mebemp.legc['LMNA',], method='spearman')
  tmp.clp.legc = our.legc[,our.annotation[colnames(our.legc)] %in% c('#6E1C71', 'blue', 'darkblue')]
  lmna.clp.mc.cors = apply(tmp.clp.legc, 1, cor, tmp.clp.legc['LMNA',], method='spearman')
  stopifnot(all(names(lmna.mebemp.cors) == names(lmna.clp.cors)))
  exported.genes = names(lmna.mebemp.cors)
  stopifnot(all(names(lmna.mebemp.mc.cors) == exported.genes))
  stopifnot(all(names(lmna.clp.mc.cors) == exported.genes))
  lmna.gene.df = data.frame(gene=exported.genes, lmna_mebemp_mc_cor=lmna.mebemp.mc.cors, lmna_clp_mc_cor=lmna.clp.mc.cors, 
                            lmna_mebemp_indiv_cor=lmna.mebemp.cors, lmna_clp_indiv_cor=lmna.clp.cors, is_in_signature = exported.genes %in% lmna.mod.fil)
  write.csv(lmna.gene.df, file=file.path(SUPP.TABLE.DIR, 's9_lmna_gene_cors.csv'), quote=F, row.names=F)

  # indiv LMNA scores
  lmna.clp.scores.unfil = colMeans(cur.exp.pooled.traj.norm.fixed.clp[lmna.mod.fil,])
  lmna.mebemp.scores.unfil = colMeans(cur.exp.pooled.traj.norm.fixed.mebemp[lmna.mod.fil,])
  tmp.exported.samps = unique(c(names(lmna.mebemp.scores.unfil), names(lmna.clp.scores.unfil)))
  tmp.exported.indivs = samp.to.indiv[tmp.exported.samps]
  lmna.scores.df = data.frame(indiv_id=tmp.exported.indivs, lmna_mebemp_score=lmna.mebemp.scores.unfil[tmp.exported.samps], lmna_clp_score=lmna.clp.scores.unfil[tmp.exported.samps])
  write.csv(lmna.scores.df, file=file.path(SUPP.TABLE.DIR, 'lmna_scores.csv'), quote=F, row.names=F)



  # and now differential expression
  indiv.id.map = get.indiv.id.map()
  indiv.to.vaf = get.indiv.to.vaf(indiv.to.used.samp.fil, indiv.id.map)
  cbc.new = get.clinical.values.new(indiv.id.map)
  for (i in 1:2) {
    cur.exp.pooled.traj.norm.fixed = all.exp.pooled.traj.norm.fixed[[i]]

    cur.samps = colnames(cur.exp.pooled.traj.norm.fixed)
    cur.indivs = samp.to.indiv[cur.samps]
    cur.samps.legc = all.samps.legcs[[i]]

    cur.male.samps = cur.samps[is.male[cur.indivs]]
    cur.female.samps = cur.samps[!is.male[cur.indivs]]
    cur.males = samp.to.indiv[cur.male.samps]
    cur.females = samp.to.indiv[cur.female.samps]

    min.exp.threshold = c(-15, -14.5)[i]
    exp.genes = rownames(cur.samps.legc)[apply(cur.samps.legc[,cur.samps], 1, max) > min.exp.threshold]
    genes.for.de = setdiff(exp.genes, (sort(names(cur.kruskal.pvals)[cur.kruskal.pvals < 1e-5])))


    tmp.cor.tests = lapply(genes.for.de, function(cur.gene) cor.test(indiv.to.vaf[cur.indivs], cur.exp.pooled.traj.norm.fixed[cur.gene, cur.samps], method='spearman'))
    tmp.wilcox.tests = lapply(genes.for.de, function(cur.gene) kruskal.test(cur.exp.pooled.traj.norm.fixed[cur.gene, cur.samps], g=indiv.to.vaf[cur.indivs] > 0))
    # both are not significant
    #tmp.cor.tests = lapply(genes.for.de, function(cur.gene) cor.test(indiv.to.vaf.tet2[cur.indivs], cur.exp.pooled.traj.norm.fixed[cur.gene, cur.samps], method='spearman'))
    #tmp.cor.tests = lapply(genes.for.de, function(cur.gene) cor.test(indiv.to.vaf.dnmt3a[cur.indivs], cur.exp.pooled.traj.norm.fixed[cur.gene, cur.samps], method='spearman'))
    tmp.cor.pvals = sapply(tmp.cor.tests, function(x) x$p.val)
    tmp.wilcox.pvals = sapply(tmp.wilcox.tests, function(x) x$p.val)
    names(tmp.cor.pvals) = genes.for.de
    names(tmp.wilcox.pvals) = genes.for.de
    head(sort(p.adjust(tmp.cor.pvals, 'BH')), n=40)
    head(sort(p.adjust(tmp.wilcox.pvals, 'BH')), n=40)

    mebemp.de.df = data.frame(gene=genes.for.de, vaf_cor=sapply(tmp.cor.tests, function(x) x$estimate), vaf_cor_pval=tmp.cor.pvals, vaf_cor_qval=p.adjust(tmp.cor.pvals, 'BH'),
                                             vaf_mann_whitney_pval=tmp.wilcox.pvals, vaf_mann_whitney_qval=p.adjust(tmp.wilcox.pvals, 'BH'))
    rownames(mebemp.de.df) = genes.for.de

    tmp.cor.tests = lapply(genes.for.de, function(cur.gene) cor.test(ages[cur.indivs], cur.exp.pooled.traj.norm.fixed[cur.gene, cur.samps], method='spearman'))
    tmp.pvals = sapply(tmp.cor.tests, function(x) x$p.val)
    names(tmp.pvals) = genes.for.de
    head(sort(p.adjust(tmp.pvals, 'BH')), n=40)

    tmp.male.cor.tests = lapply(genes.for.de, function(cur.gene) cor.test(ages[cur.males], cur.exp.pooled.traj.norm.fixed[cur.gene, cur.male.samps], method='spearman'))
    tmp.female.cor.tests = lapply(genes.for.de, function(cur.gene) cor.test(ages[cur.females], cur.exp.pooled.traj.norm.fixed[cur.gene, cur.female.samps], method='spearman'))
    tmp.male.pvals = sapply(tmp.male.cor.tests, function(x) x$p.val)
    names(tmp.male.pvals) = genes.for.de
    tmp.female.pvals = sapply(tmp.female.cor.tests, function(x) x$p.val)
    names(tmp.female.pvals) = genes.for.de
    head(sort(p.adjust(tmp.male.pvals, 'BH')), n=40)
    head(sort(p.adjust(tmp.female.pvals, 'BH')), n=40)

    mebemp.de.df$age_male_cor = sapply(tmp.male.cor.tests, function(x) x$estimate)
    mebemp.de.df$age_male_pval = tmp.male.pvals
    mebemp.de.df$age_male_qval = p.adjust(tmp.male.pvals, 'BH')
    mebemp.de.df$age_female_cor = sapply(tmp.female.cor.tests, function(x) x$estimate)
    mebemp.de.df$age_female_pval = tmp.female.pvals
    mebemp.de.df$age_female_qval = p.adjust(tmp.female.pvals, 'BH')

    tmp.wilcox.tests = lapply(genes.for.de, function(cur.gene) kruskal.test(cur.exp.pooled.traj.norm.fixed[cur.gene, cur.samps], g=is.male[cur.indivs]))
    tmp.wilcox.pvals = sapply(tmp.wilcox.tests, function(x) x$p.val)
    names(tmp.wilcox.pvals) = genes.for.de
    mebemp.de.df$sex_pvals = tmp.wilcox.pvals
    mebemp.de.df$sex_qvals = p.adjust(tmp.wilcox.pvals, 'BH')

    # male cbcs
    all.sig.genes = c()
    sig.genes.list = list()
    for (j in 1:20) {
      tmp.cor.tests = lapply(genes.for.de, function(cur.gene) cor.test(cbc.new[cur.males, j], cur.exp.pooled.traj.norm.fixed[cur.gene, cur.male.samps], method='spearman'))
      #tmp.cor.tests = lapply(genes.for.de, function(cur.gene) cor.test(cbc[cur.females, i], cur.exp.pooled.traj.norm.fixed[cur.gene, cur.females], method='spearman'))
      tmp.pvals = sapply(tmp.cor.tests, function(x) x$p.val)
      names(tmp.pvals) = genes.for.de
      print(colnames(cbc.new)[j])
      print(sum(p.adjust(tmp.pvals, 'BH') < 0.1))
      cur.sig.genes = names(tmp.pvals)[p.adjust(tmp.pvals, 'BH') < 0.1]
      all.sig.genes = c(all.sig.genes, cur.sig.genes)
      sig.genes.list[[j]] = cur.sig.genes
      mebemp.de.df[,tolower(paste0(colnames(cbc.new)[j], '_male_cor'))] = sapply(tmp.cor.tests, function(x) x$estimate)
      mebemp.de.df[,tolower(paste0(colnames(cbc.new)[j], '_male_pval'))] = tmp.pvals
      mebemp.de.df[,tolower(paste0(colnames(cbc.new)[j], '_male_qval'))] = p.adjust(tmp.pvals, 'BH')
      #(sort(p.adjust(tmp.pvals, 'BH')), n=40)
    }
    all.sig.genes = unique(all.sig.genes)
    cbc.cors = cor(cbc.new[cur.males,], t(cur.exp.pooled.traj.norm.fixed[all.sig.genes, cur.male.samps]), use='pairwise')
    shades = colorRampPalette(c("darkblue", "blue", "white", "white", "red", "darkred"))(100)
    pheatmap(cbc.cors, breaks=seq(-1, 1, length.out=100), col=shades, filename=file.path(fig.dir, 'mebemp_male_cbc.png'), width=12, height=5)

    # female cbcs
    for (j in 1:20) {
      tmp.cor.tests = lapply(genes.for.de, function(cur.gene) cor.test(cbc.new[cur.females, j], cur.exp.pooled.traj.norm.fixed[cur.gene, cur.female.samps], method='spearman'))
      tmp.pvals = sapply(tmp.cor.tests, function(x) x$p.val)
      names(tmp.pvals) = genes.for.de
      mebemp.de.df[,tolower(paste0(colnames(cbc.new)[j], '_female_cor'))] = sapply(tmp.cor.tests, function(x) x$estimate)
      mebemp.de.df[,tolower(paste0(colnames(cbc.new)[j], '_female_pval'))] = tmp.pvals
      mebemp.de.df[,tolower(paste0(colnames(cbc.new)[j], '_female_qval'))] = p.adjust(tmp.pvals, 'BH')
    }

    # export de table
    write.csv(mebemp.de.df,  file=file.path(SUPP.TABLE.DIR, c('mebemp_de_screen.csv', 'clp_de_screen.csv')[i]),  quote=F, row.names=F)

  }


}

garvan.analysis <- function(fig.dir=BASE.FIG.DIR) {
  garvan.cdata = anndata::read_h5ad(file.path(MODEL.DIR, 'pbmc_garvan_ref_cells.h5ad'))
  garvan.mdata = anndata::read_h5ad(file.path(MODEL.DIR, 'pbmc_garvan_ref_metacells.h5ad'))
  garvan.legc = t(log2(1e-5 + garvan.mdata$X / rowSums(garvan.mdata$X)))
  rownames(garvan.legc) = garvan.cdata$var$feature_name

  # CD34, hsc by age
  mc.assign.garvan = paste0('mc', garvan.cdata$obs$metacell)
  names(mc.assign.garvan) = rownames(garvan.cdata$obs)
  cell.to.indiv = garvan.cdata$obs$individual
  names(cell.to.indiv) = rownames(garvan.cdata$obs)
  garvan.cdata.obs = garvan.cdata$obs
  indiv.to.age = sapply(as.character(sort(unique(cell.to.indiv))), function(cur.indiv) {
    cur.ages = garvan.cdata.obs$age[garvan.cdata.obs$indiv == cur.indiv]
    stopifnot(all(cur.ages == cur.ages[1]))
    return(cur.ages[1])
  })

  cd34.mcs = paste0('mc', colnames(garvan.legc)[garvan.legc['CD34',] > -14.3])
  set.seed(42)
  rand.cells = unlist(lapply(unique(cell.to.indiv), function(cur.indiv) {
    candidate.cells = names(cell.to.indiv)[cell.to.indiv == cur.indiv]
    if (length(candidate.cells) < 800) {
      return(NULL)
    } else {
      return(sample(candidate.cells, 800, replace=F))
    }
  }))
  cell.to.indiv.fil = cell.to.indiv[names(cell.to.indiv) %in% rand.cells]
  cd34.frac.per.decade.test = tapply(names(cell.to.indiv.fil), floor(indiv.to.age[cell.to.indiv.fil] / 10) * 10, function(cur.cells) {
    print('decade:')
    print(unique(floor(indiv.to.age[cell.to.indiv.fil[cur.cells]] / 10) * 10))
    print('number of cells:')
    print(length(cur.cells))
    print('number of indivs:')
    print(length(unique(cell.to.indiv.fil[cur.cells])))
    #binom.test(sum(mc.assign.garvan[cur.cells] %in% cd34.mcs), sum(!(mc.assign.garvan[cur.cells] %in% c('mc-1', 'mc-2'))))
    binom.test(sum(mc.assign.garvan[cur.cells] %in% cd34.mcs), length(cur.cells))
  })
  cd34.estimates = sapply(cd34.frac.per.decade.test, function(x) x$estimate)
  cd34.lower.conf = sapply(cd34.frac.per.decade.test, function(x) x$conf.int[1])
  cd34.upper.conf = sapply(cd34.frac.per.decade.test, function(x) x$conf.int[2])
  cd34.ylim = c(0, max(cd34.estimates[1:9], cd34.lower.conf[1:9], cd34.upper.conf[1:9]))

  # fig 3D
  png(file.path(fig.dir, 'garvan_cd34_age.png'), height=800, width=400)
  log2.lim = range(log2(c(cd34.estimates, cd34.lower.conf, cd34.upper.conf)))
  plot(1:9, log2(cd34.estimates[1:9]), pch=19, ylim=log2.lim, cex=3)
  segments(1:9, log2(cd34.lower.conf[1:9]), 1:9, log2(cd34.upper.conf[1:9]), lwd=3)
  dev.off()

  # table
  garvan.indivs = unique(cell.to.indiv)
  cd34.cells = names(mc.assign.garvan)[mc.assign.garvan %in% cd34.mcs]
  sampled.cd34.cells = intersect(cd34.cells, names(cell.to.indiv.fil))
  exported.df = data.frame(indiv_id=garvan.indivs, age=indiv.to.age[garvan.indivs], decade=floor(indiv.to.age[garvan.indivs] / 10) * 10, 
             num_cells=as.character(table(cell.to.indiv)[garvan.indivs]), 
             num_cd34_pos_cells=as.character(table(cell.to.indiv[cd34.cells])[garvan.indivs]), 
	     is_selected=garvan.indivs %in% unique(cell.to.indiv.fil), num_sampled_cells=NA, num_sampled_cd34_pos_cells=NA)
  exported.df[exported.df$is_selected, 'num_sampled_cells'] = 800
  exported.df[exported.df$is_selected, 'num_sampled_cd34_pos_cells'] = as.character(table(cell.to.indiv.fil[sampled.cd34.cells])[rownames(exported.df)[exported.df$is_selected]])
  num.cd34.per.decade = tapply(as.numeric(exported.df[exported.df$is_selected, 'num_sampled_cd34_pos_cells']), exported.df[exported.df$is_selected, 'decade'], sum)
  num.cells.per.decade = tapply(as.numeric(exported.df[exported.df$is_selected, 'num_sampled_cells']), exported.df[exported.df$is_selected, 'decade'], sum)
  stopifnot(all(log2(cd34.estimates[1:9]) == log2(num.cd34.per.decade / num.cells.per.decade)))
  write.csv(exported.df, file=file.path(SUPP.TABLE.DIR, 'sx_1000_indivs_cd34.csv'), quote=F, row.names=F)
}

create.supp.tables <- function(our.cdata, our.annotation, our.assign, indiv.to.used.samp) {
  supp.tbl.dir = SUPP.TABLE.DIR
  our.cdata = anndata::read_h5ad(file.path(MODEL.DIR, '148_indiv_ref_cells.h5ad'))

  # distribution across metacells
  cell.to.meas = sprintf('%s-%s', our.cdata$obs$exp_name, our.cdata$obs$indiv_id)
  names(cell.to.meas) = rownames(our.cdata$obs)
  cell.to.samp = gsub('_\\d_ultima-', '_ultima-', gsub('_\\d-', '-', cell.to.meas))
  tmp.tbl = table(cell.to.samp, our.annotation[as.character(our.assign[names(cell.to.samp)])])
  all.samp = unique(cell.to.samp)
  indiv.to.used.samp = get.one.samp.per.indiv(all.samp, our.cdata, cell.to.samp)
  samp.to.indiv = sapply(strsplit(all.samp, '-'), function(x) x[2])
  names(samp.to.indiv) = all.samp

  age.and.sex = get.age.and.sex(our.cdata, indiv.to.used.samp, cell.to.samp)
  ages = age.and.sex$ages
  sexes = age.and.sex$sexes
  is.male = sexes == 'male'

  all.indivs = names(indiv.to.used.samp)

  indiv.id.map = get.indiv.id.map()
  tmp.cbc = get.clinical.values.new(indiv.id.map)
  tmp.cbc['N170',] = NA
  tmp.cbc['N307',] = NA
  stopifnot(all(all.indivs %in% rownames(tmp.cbc)))

  exported.indiv.md = data.frame(indiv_id=all.indivs, age=ages[all.indivs], sex=ifelse(is.male[all.indivs], 'Male', 'Female'))
  exported.cbc = tmp.cbc[all.indivs,]
  exported.cbc = cbind(data.frame(indiv_id=all.indivs), exported.cbc)

  all.muts.df = get.mutations.new(indiv.to.used.samp, indiv.id.map)
  exported.mut = filter(all.muts.df, sample %in% all.indivs)
  exported.mut$indiv_id = exported.mut$sample
  exported.mut = exported.mut[,c('indiv_id', 'Gene.refGene', 'Chr', 'Pos', 'Ref', 'Alt', 'Avg_VAF')]

  write.csv(exported.indiv.md,  file=file.path(supp.tbl.dir, 's1_individual_metadata.csv'),  quote=F, row.names=F)
  write.csv(exported.cbc,       file=file.path(supp.tbl.dir, 's2_individual_cbc.csv'),       quote=F, row.names=F)
  write.csv(exported.mut,       file=file.path(supp.tbl.dir, 's3_individual_mutations.csv'), quote=F, row.names=F)

  sync.df = read.csv(file.path(supp.tbl.dir, 'sync_scores.csv'))
  num.cells.df = read.csv(file.path(supp.tbl.dir, 'num_cells_per_ctype.csv'))
  lmna.scores.df = read.csv(file.path(supp.tbl.dir, 'lmna_scores.csv'))
  rownames(sync.df) = sync.df$indiv_id
  rownames(num.cells.df) = num.cells.df$indiv_id
  rownames(lmna.scores.df) = lmna.scores.df$indiv_id
  stopifnot(colnames(sync.df)[1] == 'indiv_id')
  stopifnot(colnames(num.cells.df)[1] == 'indiv_id')
  stopifnot(colnames(lmna.scores.df)[1] == 'indiv_id')
  indiv.scores.df = cbind(num.cells.df[all.indivs,], lmna.scores.df[all.indivs, -1], sync.df[all.indivs, -1])
  colnames(indiv.scores.df)[ncol(indiv.scores.df)] = 'sync_score'
  write.csv(indiv.scores.df, file=file.path(supp.tbl.dir, 's11_individual_scores.csv'), quote=F, row.names=F)
}

