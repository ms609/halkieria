library('TreeSearch')
# Load data from locally downloaded copy of MorphoBank matrix
nexusFile <- MorphoBank::MostRecentNexus()
my_data <- ReadAsPhyDat(nexusFile)
ignored_taxa <- c('Conotheca', 'Maxilites', 'Pauxillites',
                  'Alisina', 'Glyptoria', 'Nisusia_sulcata',
                  'Kutorgina_chengjiangensis', 'Tomteluva_perturbata',
                  'Salanygolina', 'Coolinia_pecten', 'Antigonambonites_planus',
                  'Askepasma_toddense', 'Siphonobolus_priscus',
                  'Acanthotretella_spinosa', 'Clupeafumosus_socialis',
                  'Pelagodiscus_atlanticus', 'Botsfordia', 'Eoobolus',
                  'Ussunia', 'Craniops', 'Paramicrocornus', 'Bactrotheca'
)
my_data[ignored_taxa] <- NULL
iw_data <- PrepareDataIW(my_data)
outgroup <- 'Yilingia_spiciformis'

nj.tree <- NJTree(my_data)
rooted.tree <- EnforceOutgroup(nj.tree, outgroup)
start.tree <- TreeSearch(tree=rooted.tree, dataset=my_data, maxIter=3000,
                         EdgeSwapper=RootedNNISwap, verbosity=0)

ew.tree <- Ratchet(start.tree, my_data, verbosity=3L,
                   ratchHits = 20, ratchIter = 4000,
                   searchHits = 55, searchIter = 6000,
                   swappers=list(RootedTBRSwap, RootedSPRSwap, RootedNNISwap))
ew.consensus <- RatchetConsensus(ew.tree, my_data, nSearch = 250,
                                 searchHits = 85,
                                 swappers=list(RootedTBRSwap, RootedNNISwap),
                                 verbosity=2L)

ape::write.nexus(ew.consensus, file=paste0(collapse='', "TreeSearch/hk_ew_",
                                      Fitch(ew.tree, my_data), ".nex"))
