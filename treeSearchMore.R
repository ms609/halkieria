k <- 24
library('TreeSearch')

kTreePattern <- paste0('hk_iw_k', gsub('\\.', '\\\\.', k),
                       '_\\d+\\.?\\d*\\.all\\.nex')
kFile <- list.files('TreeSearch', pattern=kTreePattern, full.names=TRUE)
bestYet <- ape::read.nexus(kFile[which.max(vapply(kFile, ApeTime, double(1)))])
if (class(bestYet) == 'multiPhylo') {
  bestYet <- bestYet[[1]]
}

# Load data from locally downloaded copy of MorphoBank matrix
nexusFile <- MorphoBank::MostRecentNexus()
myData <- ReadAsPhyDat(nexusFile)
ignoredTaxa <- c('Conotheca', 'Maxilites', 'Pauxillites',
                  'Alisina', 'Glyptoria', 'Nisusia_sulcata',
                  'Kutorgina_chengjiangensis', 'Tomteluva_perturbata',
                  'Salanygolina', 'Coolinia_pecten', 'Antigonambonites_planus',
                  'Askepasma_toddense', 'Siphonobolus_priscus',
                  'Acanthotretella_spinosa', 'Clupeafumosus_socialis',
                  'Pelagodiscus_atlanticus', 'Botsfordia', 'Eoobolus',
                  'Ussunia', 'Craniops', 'Paramicrocornus', 'Bactrotheca'
)
myData[ignoredTaxa] <- NULL
iwData <- PrepareDataIW(myData)
bestYet <- ape::drop.tip(bestYet, ignoredTaxa)

iwTree <- IWRatchet(bestYet, iwData, concavity=k,
                     ratchHits = 2, ratchIter = 4000,
                     searchHits = 32, searchIter = 4800,
                     swappers=list(RootedTBRSwap, RootedSPRSwap, RootedNNISwap),
                     verbosity=4L)
score <- IWScore(iwTree, iwData, concavity=k)
# Write a single best tree
ape::write.nexus(iwTree,
                 file=paste0("TreeSearch/hk_iw_k", k, "_",
                             signif(score, 5), ".nex", collapse=''))

iwConsensus <- IWRatchetConsensus(iwTree, iwData, concavity = k,
                                   swappers = list(RootedTBRSwap, RootedNNISwap),
                                   searchHits = 55, searchIter = 6000,
                                   nSearch = 125, verbosity = 4L)
ape::write.nexus(
  iwConsensus,
  file = paste0("TreeSearch/hk_iw_k", k, "_",
                signif(IWScore(iwConsensus[[1]], iwData, concavity=k), 5),
                ".all.nex", collapse='')
)
