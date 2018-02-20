library(MotifDb)
library(MotIV)
library(seqLogo)
length(MotifDb)
sort(table(values(MotifDb)$dataSource),decreasing=TRUE)
sort(table(values(MotifDb)$organism),decreasing=TRUE)

# http://www.nature.com/nrg/journal/v5/n4/full/nrg1315.html
# https://davetang.org/muse/2013/10/01/position-weight-matrix/
(query(query(MotifDb,'protein binding microarray'),'Scerevisiae'))
