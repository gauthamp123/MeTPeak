.metpeak.process <- function(IP,INPUT,batch_id,minimal_counts_in_fdr=10){
   print("Entered .metpeak.process()")
#   using approximate Newton's method to call peaks per gene (the fastest method)
  Ng = unique(batch_id)
  nip <- ncol(IP)
  nin <- ncol(INPUT)
  INPUT_mean <- rowMeans(INPUT)
  IP_mean <- rowMeans(IP)
  print("Step 1: Parsing input...")
  if (nip > nin) {
    avg_input <- round(matrix(rep(INPUT_mean,nip-nin),ncol=nip-nin))
    INPUT <- cbind(INPUT,avg_input) 
  }
  else if (nip < nin){
    avg_ip <- matrix(rep(IP_mean,nin-nip),ncol=nin-nip)
    IP <- cbind(IP,avg_ip) 
  }
  # initialize the variables
  pvalues <- rep(1,nrow(IP))
  print("🔄 Step 2: Looping over batches or windows...")
  for (ii in Ng){
    
    # print(ii)
    flag <- batch_id==ii
    ip=as.matrix(IP[flag,])
    input=as.matrix(INPUT[flag,])    
    # to avoid 0 read count in the input     
    bA <- max(apply(input,2,median)) + 1
    
    # to avoid ip reads is smaller than a bin of input
    # using approximate Newton written in R to do the peak calling
      res <-.betabinomial.hmm(ip,input+bA)
      cl <- .betabinomial.vit(ip,input+bA,res[[2]],res[[1]])
      clust_mean = res[[1]][1,]/colSums(res[[1]])
      maxID = which.max(clust_mean)
      
      # Peak region should have mean ratio over 0.5, meaning IP reads should be more than Input reads in peak region
      if (clust_mean[maxID] >= 0.4){
        peak_loci <- which(cl$class==maxID)
        peak_loci_end <- length(peak_loci) # no peak condition included
        if (peak_loci_end > 0){
          peak_id <- c(0,which( (peak_loci[-1] - peak_loci[-peak_loci_end] ) > 1 ),peak_loci_end)
          for (jj in 1:(length(peak_id)-1)){
            jjpeak <- peak_loci[peak_id[jj]+1]:peak_loci[peak_id[jj+1]]
            res[[3]][jjpeak,maxID] <- median(res[[3]][jjpeak,maxID])  #res[[3]] = res$postprob
          }
        }  
        pvalues[flag] = 1 - res[[3]][,maxID]
      }
      else{
        pvalues[flag] = 1
      }
  }
  
  log.fdr=log(p.adjust(pvalues,method='fdr'))
  

  ID=which( (IP_mean+INPUT_mean) > minimal_counts_in_fdr) #should be vector not matrix
  log.fdr_sig=log(p.adjust(pvalues[ID], method = "fdr"))
  log.fdr[ID]=log.fdr_sig

  log.fc=log(IP_mean/(INPUT_mean+1))

  # output result
  PW=list(log.p=log(pvalues),log.fdr=log.fdr,log.fc=log.fc)
  print("Step N: Returning to .peak.call.module()")
  return(PW)
  
}
