# This function should be replaced by the equivalent in TreeSearch 0.1.3.
GetSplits <- function (trees, tipIndex) {
  nTip <- length(tipIndex)
  if (class(trees) == 'phylo') trees <- list(trees)
  table(vapply(trees, function (tr) {
    vapply(Descendants(tr, nTip + seq_len(nTip - 1L), type='tips'),
           SplitNumber, character(2), tr, tipIndex)
  }, character(2 * (nTip - 1L))))
}

as.multiPhylo <- phytools::as.multiPhylo

GetJacks <- function (jackFile) {
  jackLines <- readLines(jackFile)
  jackTree <- TNTText2Tree(jackLines[3])
  jackTipOrder <- order(as.integer(jackTree$tip.label) + 1L)
  jackNodeOrder <- unique(unlist(Ancestors(jackTree, jackTipOrder)))[-1]
  nTntNode <- jackTree$Nnode

  treeFile <- gsub("\\.sym$", ".tre", jackFile)
  tntTrees <- ReadTntTree(treeFile, relativePath='.')
  tipLabel <- if (class(tntTrees) == 'phylo') tntTrees$tip.label else tntTrees[[1]]$tip.label
  jackTree$tip.label <- tipLabel

  jackScores <- trimws(gsub("ttag \\+\\d+ (.*); *", "\\1",
                            jackLines[length(jackLines) - (nTntNode - 2L):0]))[order(jackNodeOrder)]
  return (list(freq=gsub("^(\\d+)/.*", "\\1", jackScores),
              gc=gsub("^\\d+/(\\[?\\d+\\]?)$", "\\1", jackScores),
              tree=jackTree))
}
