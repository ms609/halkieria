#pal <- TreeSearch::brewer[[11]]
pal <- c("#F07894", "#E18758", "#C49800", "#99A700", "#55B23A", "#00B97B", "#00BAAE", "#00B3D7", "#52A2EF", "#AF8BF2", "#DE77DD") # <- colorspace::rainbow_hcl(11, 80, 65, 0, 308)

# Order matters; it'll be used (in reverse) in textFigure.R
groupCol <- c(
  Problematica = '#aaaaaa',
  Orthothecides = '#777777',
  Hyolithides = '#777777',
  Lingulids = pal[8],
  Discinids = pal[8],
  Chileids = pal[9],
  Kutorginiids = pal[9],
  Rhynchonellates = pal[9],
  Craniiforms = pal[10],
  Paterinids = pal[10],

  Aculifera = pal[6],
  Conchifera = pal[5],
  Phoronids = pal[2],
  Tommotiids = pal[3],
  Outgroup = pal[11]
)

taxonColour <- c(
  Namacalathus                    = as.character(groupCol['Problematica']),
  Yilingia                        = as.character(groupCol['Problematica']),
  Yilingia_spiciformis            = as.character(groupCol['Problematica']),
  Flustra                         = as.character(groupCol['Outgroup']),
  Amathia                         = as.character(groupCol['Outgroup']),
  Loxosomella                     = as.character(groupCol['Outgroup']),
  Serpula                         = as.character(groupCol['Outgroup']),
  Capitella                       = as.character(groupCol['Outgroup']),
  Canadia_spinosa                 = as.character(groupCol['Outgroup']),
  Canadia                         = as.character(groupCol['Outgroup']),
  Sipunculus                      = as.character(groupCol['Outgroup']),
  Dentalium                       = as.character(groupCol['Conchifera']),
  Mytilus                         = as.character(groupCol['Conchifera']),
  Pelagiella                      = as.character(groupCol['Conchifera']),
  Haliotis                        = as.character(groupCol['Conchifera']),
  Conocardium_elongatum           = as.character(groupCol['Conchifera']),
  Conocardium                     = as.character(groupCol['Conchifera']),
  Pojetaia                        = as.character(groupCol['Conchifera']),
  Pojetaia_runnegari              = as.character(groupCol['Conchifera']),
  Neopilina                       = as.character(groupCol['Conchifera']),
  Mopalia                         = as.character(groupCol['Aculifera']),
  Leptochiton                     = as.character(groupCol['Aculifera']),
  Tonicella                       = as.character(groupCol['Aculifera']),
  Wirenia                         = as.character(groupCol['Aculifera']),
  Chaetoderma                     = as.character(groupCol['Aculifera']),
  Polysacos                       = as.character(groupCol['Aculifera']),
  Polysacos_vickersianum          = as.character(groupCol['Aculifera']),
  Glaphurochiton                  = as.character(groupCol['Aculifera']),
  Glaphurochiton_carbonarius      = as.character(groupCol['Aculifera']),
  Calvapilosa                     = as.character(groupCol['Aculifera']),
  Calvapilosa_kroegeri            = as.character(groupCol['Aculifera']),
  Acaenoplax_hayae                = as.character(groupCol['Aculifera']),
  Acaenoplax                      = as.character(groupCol['Aculifera']),
  Kulindroplax_perissokomos       = as.character(groupCol['Aculifera']),
  Kulindroplax                    = as.character(groupCol['Aculifera']),
  Phthipodochiton_thraivensis     = as.character(groupCol['Aculifera']),
  Phthipodochiton                 = as.character(groupCol['Aculifera']),
#  Halkieria                       = as.character(groupCol['Problematica']),
  Halkieria                       = as.character('#000000'),
  Halkieria_evangelista           = 'black',
#  Halkieria_evangelista           = as.character(groupCol['Problematica']),
#  Wiwaxia_corrugata               = as.character(groupCol['Problematica']),
#  Wiwaxia                         = as.character(groupCol['Problematica']),
  Wiwaxia_corrugata               = 'black',
  Wiwaxia                         = 'black',
  Orthrozanclus                   = as.character(groupCol['Problematica']),
  Odontogriphus                   = as.character(groupCol['Problematica']),
  Odontogriphus_omalus            = as.character(groupCol['Problematica']),
  Cotyledion_tylodes              = as.character(groupCol['Problematica']),
  Cotyledion                      = as.character(groupCol['Outgroup']),
  Novocrania                      = as.character(groupCol['Craniiforms']),
  Craniops                        = as.character(groupCol['Craniiforms']),
  Ussunia                         = as.character(groupCol['Craniiforms']),
  Gasconsia                       = as.character(groupCol['Craniiforms']),
  Heliomedusa_orienta             = as.character(groupCol['Tommotiids']),
  Heliomedusa                     = as.character(groupCol['Tommotiids']),
  Micrina                         = as.character(groupCol['Tommotiids']),
  Mickwitzia_muralensis           = as.character(groupCol['Tommotiids']),
  Mickwitzia                      = as.character(groupCol['Tommotiids']),
  Micromitra                      = as.character(groupCol['Paterinids']),
  Askepasma_toddense              = as.character(groupCol['Paterinids']),
  Askepasma                       = as.character(groupCol['Paterinids']),
  Pelagodiscus_atlanticus         = as.character(groupCol['Discinids']),
  Pelagodiscus                    = as.character(groupCol['Discinids']),
  Mummpikia_nuda                  = as.character(groupCol['Problematica']),
  Mummpikia                       = as.character(groupCol['Problematica']),
  Lingula                         = as.character(groupCol['Lingulids']),
  Eoobolus                        = as.character(groupCol['Lingulids']),
  Botsfordia                      = as.character(groupCol['Lingulids']),
  Siphonobolus_priscus            = as.character(groupCol['Lingulids']),
  Siphonobolus                    = as.character(groupCol['Lingulids']),
  Clupeafumosus_socialis          = as.character(groupCol['Lingulids']),
  Clupeafumosus                   = as.character(groupCol['Lingulids']),
  Lingulellotreta_malongensis     = as.character(groupCol['Lingulids']),
  Lingulellotreta                 = as.character(groupCol['Lingulids']),
  Acanthotretella_spinosa         = as.character(groupCol['Lingulids']),
  Acanthotretella                 = as.character(groupCol['Lingulids']),
  Lingulosacculus                 = as.character(groupCol['Problematica']),
  Phoronis                        = as.character(groupCol['Phoronids']),
  Dailyatia                       = as.character(groupCol['Tommotiids']),
  Eccentrotheca                   = as.character(groupCol['Tommotiids']),
  Yuganotheca_elegans             = as.character(groupCol['Problematica']),
  Yuganotheca                     = as.character(groupCol['Problematica']),
  Longtancunella_chengjiangensis  = as.character(groupCol['Problematica']),
  Longtancunella                  = as.character(groupCol['Problematica']),
  Paterimitra                     = as.character(groupCol['Tommotiids']),
  Paramicrocornus                 = as.character(groupCol['Hyolithides']),
  Pedunculotheca_diania           = as.character(groupCol['Hyolithides']),
  Pedunculotheca                  = as.character(groupCol['Hyolithides']),
  Haplophrentis_carinatus         = as.character(groupCol['Hyolithides']),
  Haplophrentis                   = as.character(groupCol['Hyolithides']),
  Maxilites                       = as.character(groupCol['Hyolithides']),
  Pauxillites                     = as.character(groupCol['Hyolithides']),
  Bactrotheca                     = as.character(groupCol['Orthothecides']),
  Conotheca                       = as.character(groupCol['Orthothecides']),
  Cupitheca_holocyclata           = as.character(groupCol['Orthothecides']),
  Cupitheca                       = as.character(groupCol['Orthothecides']),
  Tomteluva_perturbata            = as.character(groupCol['Problematica']),
  Tomteluva                       = as.character(groupCol['Problematica']),
  Salanygolina                    = as.character(groupCol['Chileids']),
  Coolinia_pecten                 = as.character(groupCol['Chileids']),
  Coolinia                        = as.character(groupCol['Chileids']),
  Antigonambonites_planus         = as.character(groupCol['Chileids']),
  Antigonambonites                = as.character(groupCol['Chileids']),
  Kutorgina_chengjiangensis       = as.character(groupCol['Kutorginiids']),
  Kutorgina                       = as.character(groupCol['Kutorginiids']),
  Nisusia_sulcata                 = as.character(groupCol['Kutorginiids']),
  Nisusia                         = as.character(groupCol['Kutorginiids']),
  Alisina                         = as.character(groupCol['Problematica']),
  Glyptoria                       = as.character(groupCol['Rhynchonellates']),
  Orthis                          = as.character(groupCol['Rhynchonellates']),
  Terebratulina                   = as.character(groupCol['Rhynchonellates']))


ColMissing <- function (omit, position='bottomleft') {
  if (length(omit) > 0)
    MarkMissing(omit, position, text.font=3, cex=0.8, text.col=taxonColour[omit])
}

WriteNumber <- function (n) {
  if (n == 0) {
    'zero'
  } else if (0 < n && n < 10) {
    c('one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine')[n]
  } else {
    n
  }
}

UnitEdges <- function (tree) {
  tree$edge.length <- rep(1, dim(tree$edge)[1])
  # Return:
  tree
}

SetPar <- function() par(mar=rep(0.2, 4), cex=0.8)

ColPlot <- function (tree) {
  SetPar()
  plot(tree, tip.color = taxonColour[tree$tip.label], lwd=2)
}

# Plot results for each value of k
PlotPanel <- function (treeList, i, jackFile = NULL) {
  if (i > length(treeList) || length(treeList[[i]]) == 0) {
    cat("\n > Results not available for panel", i)
  } else {
    ColPlot(RootTree(consensus(treeList[[i]]), rootingTips))
    text(1, 1, paste0('k = ', kValues[i]), pos=4)
    if (!is.null(jackFile)) {
      jacks <- GetJacks(jackFile)
      nodeScores <- as.integer(jacks$freq)
      nodelabels(paste0(c('', nodeScores), "\n"),
                 col=c('black', SupportColour(nodeScores / 200L + 0.5, show1=FALSE)),
                 adj=0, pos=2, frame='none', cex=0.75)
    }
  }
}

LabelNodes <- function (x, col=1, line=1) {
  nodelabels(paste(rep("\n\n", 2 - line), x, rep("\n\n", line)),
             frame='none', pos=2, cex=0.8, col=brewer[[4]][col])
}
