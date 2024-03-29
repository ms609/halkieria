library('TreeSearch')
# Load data from locally downloaded copy of MorphoBank matrix
nexusFile <- MorphoBank::MostRecentNexus()
my_data <- ReadAsPhyDat(nexusFile)
ignoredTaxa <- c('Conotheca', 'Maxilites', 'Pauxillites', 'Alisina',
                 'Amathia', # There because it has setae; Flustra doesn't
                 'Yuganotheca_elegans', 'Lingulellotreta_malongensis',
                 'Glyptoria', 'Nisusia_sulcata', 'Kutorgina_chengjiangensis',
                 'Tomteluva_perturbata', 'Salanygolina', 'Coolinia_pecten',
                 'Antigonambonites_planus', 'Askepasma_toddense',
                 'Siphonobolus_priscus', 'Acanthotretella_spinosa',
                 'Longtancunella_chengjiangensis', 'Lingulosacculus',
                 'Mummpikia_nuda', 'Cupitheca_holocyclata',
                 'Clupeafumosus_socialis', 'Botsfordia', 'Eoobolus', 'Ussunia',
                 'Craniops', 'Paramicrocornus', 'Bactrotheca')
my_data[ignoredTaxa] <- NULL
iw_data <- PrepareDataIW(my_data)
outgroup <- 'Yilingia_spiciformis'

nj.tree <- NJTree(my_data)
rooted.tree <- EnforceOutgroup(nj.tree, outgroup)
start.tree <- TreeSearch(tree=rooted.tree, dataset=my_data, maxIter=3000,
                         EdgeSwapper=RootedNNISwap, verbosity=0)

ew.tree <- Ratchet(start.tree, my_data, verbosity=3L,
                   ratchHits = 42, ratchIter = 1000,
                   searchHits = 55, searchIter = 12000,
                   swappers=list(RootedTBRSwap, RootedSPRSwap, RootedNNISwap))
ew.consensus <- RatchetConsensus(ew.tree, my_data, nSearch = 250,
                                 searchHits = 85,
                                 swappers=list(RootedTBRSwap, RootedNNISwap),
                                 verbosity=2L)

ape::write.nexus(ew.consensus, file=paste0(collapse='', "TreeSearch/hk_ew_",
                                      Fitch(ew.tree, my_data), ".nex"))
