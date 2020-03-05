#Creates a simulated bulk data set from sc data
#
#
#
bulkData<-function(scData, nSample=30, nSC = 25) {
  bulk_exp = matrix(nrow = nSample, ncol = ncol(scData))
  colnames(bulk_exp) = colnames(scData)
  for(i in (1:nSample)) {
    bulk_exp[i,] = colMeans(scData[sample(nrow(scData),nSC), ])
    
  }
  return(bulk_exp)
}

#Pulls any in a network adjacency matrix interactions above a threshold out and returns
#a matrix of transcription factor and target gene
#
#
#
getInteractions<-function(net, threshold = 0.3) {
  index = which(net > threshold, arr.ind = TRUE)
  grn = matrix(ncol = 2, nrow = nrow(index))
   grn[,1] = colnames(net)[index[,2]]
   grn[,2] = rownames(net)[index[,1]]
  colnames(grn) = c("TF","Target")
  return(grn)
}

#Made to assess a reconstructed GRN when given a gold standard
#Turns out is obsolete because minet has a very nice function to do this for us
#
#
grnAssess<-function(grn,grn_gs, net, verbose = FALSE) {
  tp = nrow(rbind(grn,grn_gs))-nrow(unique(rbind(grn,grn_gs)))
  fp = nrow(grn) - tp
  tn = (ncol(net)*nrow(net)) - (nrow(grn_gs) - tp)
  fn = (ncol(net)*nrow(net)) - tn
  precision = tp/(tp+fp)
  recall = tp/(tp+fn)
  
  if (verbose) {
    cat('True Positives ', tp, '\n')
    cat('False Positives ', fp, '\n')
    cat('True Negatives ', tn, '\n')
    cat('False Negatives ', fn, '\n')
    cat('Precission: ',precision, '\n')
    cat('Recall: ', recall, '\n')
    return(c(precision,recall))
  }
  else{
    return(c(precision,recall))
  }
}


#Takes the full network of data and reduces the column dimension to just include a list of TFs
#
#
#
#
reduceNet<-function(net,tf_list) {
  new_net = c()
  n = 0
  for(i in (1:length(tf_list))) {
    if(tf_list[i] %in% rownames(net)) {
      n = n + 1 
      new_net = cbind(new_net,net[,which(rownames(net) == tf_list[i])])
      colnames(new_net)[n] = tf_list[i]
    }
  }
  return(new_net)
}

#Wrapper function for the above function
#Will return the reconstructed grn and output the precision adn recall curve if PR is set to TRUE
#
#
reconstruct<-function(sc_exp, TF_list, gs_grn, thresh = 0.3, normalize = TRUE, PR = FALSE) {
  pr = c()
  mim <- build.mim(sc_exp,estimator="pearson")
  net <- clr(mim)
  net = reduceNet(net,TF_list)
  n = 0
  if (normalize == TRUE) {
    net = sqrt(net)
    #capture.output(net = as.matrix(NormalizeData(net)))
    
    net = (net - min(net))/(max(net)-min(net))
  }
  if (PR == TRUE) {
    for(i in (0:99)) {
      n = n + 1
      if (i %% 10 == 0) {
        cat('...', i, ' %...',)
      }
      grn = getInteractions(net, threshold=(i/100))
      pr = rbind(pr, grnAssess(grn,gs_grn,net))
      rownames(pr)[n] = ((i/100))
    }
  }
  else{
    grn = getInteractions(net, threshold = thresh)
  }
  grnAssess(grn,gs_grn,net)
  
  colnames(pr) = c('Precision','Recall')
  plot((pr[,2]), (pr[,1]) , main = "Precision Recall Curve", ylab = "Precision", xlab = "Recall")
  View(pr)
  return(grn)
  
}
makeAdj<-function(net, gs_grn) {
  adj = matrix(data = 0,nrow=nrow(net),ncol=ncol(net))
  colnames(adj) = colnames(net)
  rownames(adj) = rownames(net)
  for(i in (1:nrow(gs_grn))) {
    index1 = which(colnames(net) == gs_grn[i,1])
    index2 = which(rownames(net) == gs_grn[i,2])
    adj[index2,index1] = 1
  }
  return(adj)
}


 