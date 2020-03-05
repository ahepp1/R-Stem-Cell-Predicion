expDat3 <- AverageExpression(bertie_exp)
expDat3 <- data.matrix(expDat3[[1]])
interactionScores <- matrix(0L,nrow = ncol(expDat3), ncol = ncol(expDat3))
isEmpty <- function(x) return(length(x) ==0 )
meanLigandExpression <- 0
ligandCount <- 0
meanReceptorExpression <- 0
receptorCount <- 0
rownames <- row.names(expDat)
for(i in 1:ncol(expDat3)) {
  for(j in 1:length(expDat3[,i])) {
    if(length(intersect(rownames[j], ligands)) == 1) {
      meanLigandExpression <- meanLigandExpression + expDat3[j,i]
      ligandCount <- ligandCount + 1
    } 
    if(length(intersect(rownames[j], receptors)) == 1) {
      meanReceptorExpression <- meanReceptorExpression + expDat3[j,i]
      receptorCount <- receptorCount + 1
    }
  }
}
meanLigandExpression <- meanLigandExpression/ligandCount
meanReceptorExpression <- meanReceptorExpression/receptorCount