cname.to.exp.name <- function(cnames) {
  return(sapply(strsplit(cnames, '_'), function(x) paste0(x[1:(length(x) - 1)], collapse='_')))
}

get.tech.replicates <- function(cell.to.meas, cnames) {
  cell.to.meas = cell.to.meas[cnames]
  all.meas = unique(cell.to.meas)

  tech.reps = do.call(rbind, lapply(all.meas, function(cur.meas) {
    splitted = strsplit(cur.meas, '-')[[1]]
    cur.exp = splitted[1]
    cur.indiv = splitted[2]
    potential.rep = sprintf('%s_ultima-%s', cur.exp, cur.indiv)
    if(potential.rep %in% all.meas) {
      return(list(cur.meas, potential.rep))
    } else {
      return(NULL)
    }
  }))
  return(tech.reps)
}


get.tech.rep.pools <- function(cdata.x, cdata.obs, cnames) {
  cell.to.meas = sprintf('%s-%s', cdata.obs$exp_name, cdata.obs$indiv_id)
  names(cell.to.meas) = rownames(cdata.obs)
  tech.reps = get.tech.replicates(cell.to.meas, cnames)

  common.cells.per.meas = lapply(1:nrow(tech.reps), function(i) {
    cur.illum.cnames = names(cell.to.meas)[cell.to.meas == tech.reps[i, 1]]
    cur.ultima.cnames = names(cell.to.meas)[cell.to.meas == tech.reps[i, 2]]
    cur.illum.barcodes = cdata.obs[cur.illum.cnames, 'barcode']
    cur.ultima.barcodes = cdata.obs[cur.ultima.cnames, 'barcode']
    cur.common.barcodes = intersect(cur.illum.barcodes, cur.ultima.barcodes)
    return(list(cur.illum.cnames[cur.illum.barcodes %in% cur.common.barcodes],
                cur.ultima.cnames[cur.ultima.barcodes %in% cur.common.barcodes]))
  })

  tmp.cnames = unlist(common.cells.per.meas)
  tmp.exp.names = cname.to.exp.name(tmp.cnames)
  tmp.pools = do.call(cbind, tapply(tmp.cnames, tmp.exp.names, function(cur.cnames) {
    cur.pool = colSums(cdata.x[cur.cnames,])
    return(cur.pool)
  }))
  rep.illum.pool = tmp.pools[,!grepl('ultima', colnames(tmp.pools))]
  rep.ultima.pool = tmp.pools[,grepl('ultima', colnames(tmp.pools))]
  return(list(rep.illum.pool=rep.illum.pool, rep.ultima.pool=rep.ultima.pool))
}


