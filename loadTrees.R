library('ape')
source('treeFunctions.R')
ew.best <- list.files('TreeSearch', pattern='hk_ew_\\d*\\.nex', full.names=TRUE)
ew.trees <- read.nexus(file=ew.best[which.max(file.mtime(ew.best))])
ew.trees <- if (class(ew.trees) == 'multiPhylo') unique(ew.trees) else ew.trees
tipIndex <- sort(ew.trees[[1]]$tip.label)
ew.trees <- as.multiPhylo(ew.trees)
ew.trees <- lapply(ew.trees, RenumberTips, tipIndex)
ew.trees <- lapply(ew.trees, SortTree)
class(ew.trees) <- 'multiPhylo'
iw.trees <- lapply(kValues, function (k) {
  iw.best <- list.files('TreeSearch',
                        pattern=paste0('hk_iw_k',
                                       gsub('\\.', '\\\\.', k),
                                       '_\\d+\\.?\\d*\\.all\\.nex'),
                        full.names=TRUE)
  # Return:
  if (length(iw.best) == 0) {
    list()
  } else {
    loadedTrees <- read.nexus(iw.best[which.max(vapply(iw.best, ApeTime, double(1)))])
    loadedTrees <- if (class(loadedTrees) == 'multiPhylo') {
      unique(loadedTrees)
    } else {
      as.multiPhylo(loadedTrees)
    }
    loadedLabels <- loadedTrees[[1]]$tip.label
    if (all(tipIndex %in% loadedLabels) && all(loadedLabels %in% tipIndex)) {
      ret <- lapply(loadedTrees, RenumberTips, tipIndex)
      ret <- structure(lapply(ret, SortTree), class= 'multiPhylo')
      ret
    } else {
      warning("Tip labels do not match ew labels for k = ", k, '.',
              '\n  ew tree has ', length(tipIndex), ' tips, iw tree ',
              length(loadedLabels))

      list()
    }
  }
})
iw.treesLoaded <- vapply(iw.trees, length, 0) > 0
iw.exist <- iw.trees[iw.treesLoaded]

allTrees <- ew.trees
for (i in seq_along(iw.exist)) {
  allTrees <- c(allTrees, iw.exist[[i]])
}
