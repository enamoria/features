# This scripts used package 'protr' to calculate SOCN and QSO

scriptPath <- function() {
  getSrcDirectory(scriptPath);
}
scriptPath()
setwd(scriptPath())

install.packages('protr')
library('protr')

# x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))[[1]]
x = readFASTA(system.file('protseq/P00750.fasta', package = 'protr'))
x
extractSOCN(x)

