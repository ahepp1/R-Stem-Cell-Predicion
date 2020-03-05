expDat <- AverageExpression(bertie_exp)
expDat <- data.matrix(expDat[[1]])
rownames <- row.names(expDat)
interactionScores2 <- matrix(0L,nrow = 11, ncol = 11)
isEmpty <- function(x) return(length(x) ==0 )
ligands <- intersect(row.names(expDat), LR_pairs$mouse.ligand)
receptors <- intersect(row.names(expDat), LR_pairs$mouse.receptor)
allLigands <- LR_pairs$mouse.ligand
allReceptors <- LR_pairs$mouse.receptor
for(k in 1:11) {
  for(i in length(expDat[,k])) {
    #if(!isEmpty(intersect(expDat[i,k], ligands)) && expDat[i,k] > 0) {
      myLigand <- intersect(expDat[i,k], ligands)
      index <- match(myLigand, allLigands)
      myReceptor <- allReceptors[index]
      receptorCells <- expDat[myReceptor,]
      receptorIndex <- match(myReceptor, rownames)
      for(j in 1:length(receptorCells)) {
        if(receptorCells[j] > 0) {
          interactionScores[i,j] = interactionScores[i,j] + expDat[i,k] * expDat[receptorIndex,j]
        }
      }
    #}
  }
}