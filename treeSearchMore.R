kValues <- k <- 10.5
library('TreeSearch')
source('loadTrees.R')

# Load data from locally downloaded copy of MorphoBank matrix
nexusFile <- MorphoBank::MostRecentNexus()
my_data <- ReadAsPhyDat(nexusFile)
ignored_taxa <- c('Conotheca', 'Maxilites', 'Pauxillites',
                  'Probactrotheca') # Also manually update tnt.run using `taxcode-`
my_data[ignored_taxa] <- NULL
iw_data <- PrepareDataIW(my_data)
outgroup <- 'Yilingia_spiciformis'

start.tree <- iw.trees[[1]][[1]]

iw.tree <- IWRatchet(start.tree, iw_data, concavity=k,
                     ratchHits = 6, ratchIter = 4000,
                     searchHits = 48, searchIter = 1600,
                     swappers=list(RootedTBRSwap, RootedSPRSwap, RootedNNISwap),
                     verbosity=4L)
score <- IWScore(iw.tree, iw_data, concavity=k)
# Write a single best tree
ape::write.nexus(iw.tree,
                 file=paste0("TreeSearch/hk_iw_k", k, "_",
                             signif(score, 5), ".nex", collapse=''))

iw.consensus <- IWRatchetConsensus(iw.tree, iw_data, concavity = k,
                                   swappers = list(RootedTBRSwap, RootedNNISwap),
                                   searchHits = 55, searchIter = 6000,
                                   nSearch = 250, verbosity = 4L)
ape::write.nexus(
  iw.consensus,
  file = paste0("TreeSearch/hk_iw_k", k, "_",
                signif(IWScore(iw.consensus[[1]], iw_data, concavity=k), 5),
                ".all.nex", collapse='')
)
