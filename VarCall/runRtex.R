library(knitr)
setwd("/research/bsi/projects/breast/s108235.tripneg_Couch/projects/s108235.tripneg/breast_requests/Predisposition_genes/Gilbert/varcall/Rtex/")

knit("mayoPreprocessing.Rtex")
cat('finsihed mayoPreprocessing.Rtex\n')
knit("nciPreprocessing.Rtex")

cat('finsihed nciPreprocessing.Rtex\n')

knit("BRCA2CombinedMave24.ldaER.Run1.Rtex")
cat('finsihed BRCA2CombinedMave25.ldaER.Run1.Rtex\n')

knit("BRCA2CombinedMave24.ldaER.Run2.Rtex")
cat('finsihed BRCA2CombinedMave24.ldaER.Run2.Rtex\n')
