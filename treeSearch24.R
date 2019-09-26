library('TreeSearch')
# Load data from locally downloaded copy of MorphoBank matrix
nexusFile <- MorphoBank::MostRecentNexus()
my_data <- ReadAsPhyDat(nexusFile)
ignored_taxa <- c('Conotheca', 'Maxilites', 'Pauxillites', 'Probactrotheca')
my_data[ignored_taxa] <- NULL
iw_data <- PrepareDataIW(my_data)
outgroup <- 'Yilingia_spiciformis'

nj.tree <- NJTree(my_data)
rooted.tree <- EnforceOutgroup(nj.tree, outgroup)
start.tree <- TreeSearch(tree=rooted.tree, dataset=my_data, maxIter = 20000,
                         EdgeSwapper=RootedNNISwap, verbosity=0)

k <- 24
iw.tree <- IWRatchet(start.tree, iw_data, concavity=k,
                     ratchHits = 20, ratchIter = 4000,
                     searchHits = 64, searchIter = 4000,
                     swappers = list(RootedTBRSwap, RootedSPRSwap, RootedNNISwap),
                     verbosity = 4L)
score <- IWScore(iw.tree, iw_data, concavity = k)
# Write a single best tree
ape::write.nexus(iw.tree,
                 file=paste0("TreeSearch/hk_iw_k", k, "_",
                             signif(score, 5), ".nex", collapse=''))

iw.consensus <- IWRatchetConsensus(iw.tree, iw_data, concavity = k,
                                   swappers = list(RootedTBRSwap, RootedNNISwap),
                                   searchHits = 55, searchIter = 8000,
                                   nSearch = 250, verbosity = 4L)
ape::write.nexus(iw.consensus,
  file=paste0("TreeSearch/hk_iw_k", k, "_",
        signif(IWScore(iw.consensus[[1]], iw_data, concavity=k), 5),
        ".all.nex", collapse=''))
