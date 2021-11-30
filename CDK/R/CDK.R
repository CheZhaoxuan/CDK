#' @title construct different kinship matrix
#'
#' @description
#' \code{G_matrix}use different methods to construct the kinship matrix
#'
#' @param geno a matrix
#' @param maf a vector of type numeric
#' @param A a matrix
#' @param EV a vector of type numeric
#' @param nPCs the number of front PCs, factor
#' @param method a character, options are VanRaden1, VanRaden2, VanRaden3, Weighted, GEMMA1, GEMMA2, GAB, Genetic Distance, GIBS, KING, HIERFSTAT, one-SNP haplotype, GAPIT1, GAPIT2
#'
#' @return a matrix
#' @export
#'
#' @examples
#' G<-CDK(geno=geno_matrix,method="VanRaden1")

CDK<-
  function(geno,maf=NULL,A=NULL,EV=NULL,nPCs=0,
           method="VanRaden1"){
    #--------------------------------------------------------------------------------------------------------#
    # Object: To calculate the G matrix in different methods												                         #
    #																										                                                     #
    # Input:	 						 																	                                                 #
    # geno: genotype in numeric format, pure 0, 1, 2 matrix; geno is n individual rows by m snps columns	   #
    # maf: a vector; maf is allele frequency of second allele                                       	       #
    # A: a matrix; A is a pedigree in numeric format                                                         #
    # EV: Eigen values                                                                                       #
    # nPCs: the number of front PCs excluded to calculate kinship                                            #			                                                                           	 #													 #
    # method: the method for calculating the G matrix;                                                      #
    #         options are VanRaden1, VanRaden2, VanRaden3, Weighted, GEMMA1, GEMMA2, GAB, Genetic Distance, GIBS, KING, HIERFSTAT, one-SNP haplotype, GAPIT1, GAPIT2						 				 #
    #--------------------------------------------------------------------------------------------------------#
    if (!is.matrix(geno)) stop("geno should be a matrix")
    if (!is.numeric(geno)) stop ("geno contains non numeric values")
    if (any(is.na(geno))) warning("geno contains NA values")
    if (is.null(maf)) maf<-colMeans(geno,na.rm=T)/2
    names(maf)<-colnames(geno)
    if (ncol(geno)!=length(maf)) stop("length of maf vector should equal number of columns in geno")
    if (sum(names(maf)!=colnames(geno))!=0) stop ("SNP names on vector maf and geno don't match")
    if (is.na(sum(maf))) stop ("allelic maf vector can't have NA values")
    if (is.null(method)) stop("method is NULL")
    if (!method %in% c("VanRaden1","VanRaden2","VanRaden3","Weighted","GEMMA1","GEMMA2","GAB","Genetic Distance","GIBS","KING","HIERFSTAT","one-SNP haplotype","GAPIT1","GAPIT2")) stop("This method is not supported")
    if ((method %in% c("VanRaden1","VanRaden2","VanRaden3","Weighted","GEMMA1","GEMMA2","GAB","Genetic Distance","GIBS","KING","HIERFSTAT","one-SNP haplotype","GAPIT1","GAPIT2")) & (length(method)!=1)) stop("Only one method can be input at a time")

    VanRaden1_matrix<-function(geno,maf){
      #--------------------------------------------------------------------------------------------------------
      # Object: To calculate the VanRaden1 matrix
      # Method source:
      #	VanRaden PM. Efficient methods to compute genomic predictions, J Dairy Sci 2008;91:4414-4423
      #
      # Input:
      # geno: genotype in numeric format, pure 0, 1, 2 matrix; geno is n individual rows by m snps columns
      # maf: a vector; maf is allele frequency of second allele
      #
      # Output:
      # VanRaden1 matrix
      #--------------------------------------------------------------------------------------------------------
      #Remove invariants
      Fa <- colMeans(geno,na.rm=T)/2
      index.re <- Fa>=1|Fa<=0
      geno <- geno[,!index.re]
      maf <- colMeans(geno,na.rm=T)/2

      #Calculate Z matrix
      M <- geno
      P <- maf
      Z <- as.matrix(M-2*P)

      #G=tcrossprod((geno), (geno))
      G <- tcrossprod(Z,Z)

      #Adjust
      Adjust <- function(x){x*(1-x)}
      Adj <- 2*sum(Adjust(P))
      VanRaden1 <- G/Adj

      return(VanRaden1)}

    VanRaden2_matrix<-function(geno){
      #--------------------------------------------------------------------------------------------------------
      # Object: To calculate the G0.5 matrix
      # Method source:
      #	Forni S, Aguilar I, Misztal I. Different genomic relationship matrices for single-step analysis using phenotypic, pedigree and genomic information, Genet Sel Evol 2011;43:1
      #
      # Input:
      # geno: genotype in numeric format, pure 0, 1, 2 matrix; geno is n individual rows by m snps columns
      #
      # Output:
      # VanRaden2 matrix
      #--------------------------------------------------------------------------------------------------------
      #Remove invariants
      Fa <- colMeans(geno,na.rm=T)/2
      index.re <- Fa>=1|Fa<=0
      geno <- geno[,!index.re]

      #Calculate Z matrix
      M <- geno
      nSNP <- ncol(geno)
      nInd <- nrow(geno)
      P <- rep(0.5,nSNP)
      Z <- as.matrix(M-2*P)

      #G=tcrossprod((geno), (geno))
      G <- tcrossprod(Z,Z)

      #Adjust
      Adj <- nSNP/2
      VanRaden2 <- G/Adj

      return(VanRaden2)}

    VanRaden3_matrix<-function(geno,A){
      #--------------------------------------------------------------------------------------------------------
      # Object: To calculate the VanRaden3 matrix
      # Method source:
      #	VanRaden PM. Efficient methods to compute genomic predictions, J Dairy Sci 2008;91:4414-4423
      #
      # Input:
      # geno: genotype in numeric format, pure 0, 1, 2 matrix; geno is n individual rows by m snps columns
      # A: a matrix; A is a pedigree in numeric format
      #
      # Output:
      # VanRaden3 matrix
      #--------------------------------------------------------------------------------------------------------
      if (is.null(A)) stop("A should be provided")
      if (!is.matrix(A))stop("A should be a matrix")
      if (!is.numeric(A)) stop ("A contains non numeric values")

      #Remove invariants
      Fa <- colMeans(geno,na.rm=T)/2
      index.re <- Fa>=1|Fa<=0
      geno <- geno[,!index.re]

      M <- geno
      nSNP <- ncol(geno)
      nInd <- nrow(geno)

      #Calculate parameter g0 and g1
      LHS11 <- nInd^2
      LHS12 <- sum(colSums(A))
      LHS21 <- LHS12
      LHS22 <- sum(colSums(A^2))
      RHS1 <- sum(colSums(tcrossprod(M,M)))
      RHS2 <- sum(colSums(tcrossprod(M,M)*A))

      LHS <- matrix(c(LHS11,LHS12,LHS21,LHS22),2,2)
      RHS <- c(RHS1,RHS2)

      g0 <- solve(LHS,RHS)[1]
      g1 <- solve(LHS,RHS)[2]

      #Calculate G matrix
      a1<-matrix(1,nInd,nSNP)
      a11<-tcrossprod(a1,a1)

      VanRaden3 <- (tcrossprod(M,M)-g0*a11)/g1

      return(VanRaden3)}

    Weighted_matrix<-function(geno,maf){
      #--------------------------------------------------------------------------------------------------------
      # Object: To calculate the Weighted matrix
      # Method source:
      #	VanRaden PM. Efficient methods to compute genomic predictions, J Dairy Sci 2008;91:4414-4423.
      #
      # Input:
      # geno: genotype in numeric format, pure 0, 1, 2 matrix; geno is n individual rows by m snps columns
      # maf: a vector; maf is allele frequency of second allele
      #
      # Output:
      # Weighted matrix
      #--------------------------------------------------------------------------------------------------------
      #Remove invariants
      Fa <- colMeans(geno,na.rm=T)/2
      index.re <- Fa>=1|Fa<=0
      geno <- geno[,!index.re]
      maf <- colMeans(geno,na.rm=T)/2

      #Calculate D matrix
      M <- geno
      m <- ncol(geno)
      P <- maf

      D_matrix <- function(x) {1/m*(2*x*(1-x))}
      D <- diag(D_matrix(P))

      #G = Z%*%D%*%t(Z)
      Z <- as.matrix(M-2*P)
      L <- Z%*%D
      Weighted <- tcrossprod(L,Z)

      return(Weighted)}

    GEMMA1_matrix<-function(geno){
      #--------------------------------------------------------------------------------------------------------
      # Object: To calculate the GEMMA1 matrix
      # Method source:
      #	Zhou X, Stephens M. Genome-wide efficient mixed-model analysis for association studies, Nat Genet 2012;44:821-824.
      #
      # Input:
      # geno: genotype in numeric format, pure 0, 1, 2 matrix; geno is n individual rows by m snps columns
      #
      # Output:
      # GEMMA1 matrix
      #--------------------------------------------------------------------------------------------------------
      #Remove invariants
      Fa <- colMeans(geno,na.rm=T)/2
      index.re <- Fa>=1|Fa<=0
      geno <- geno[,!index.re]

      #Calculate G matrix
      M <- geno
      p <- ncol(geno)

      cmean <- colMeans(M,na.rm = TRUE)
      W <- as.matrix(sweep(M,2,cmean))

      GEMMA1 <- tcrossprod(W)/p

      return(GEMMA1)}

    GEMMA2_matrix<-function(geno){
      #--------------------------------------------------------------------------------------------------------
      # Object: To calculate the GEMMA2 matrix
      # Method source:
      #	Zhou X, Stephens M. Genome-wide efficient mixed-model analysis for association studies, Nat Genet 2012;44:821-824.
      #
      # Input:
      # geno: genotype in numeric format, pure 0, 1, 2 matrix; geno is n individual rows by m snps columns
      #
      # Output:
      # GEMMA2 matrix
      #--------------------------------------------------------------------------------------------------------
      #Remove invariants
      Fa <- colMeans(geno,na.rm=T)/2
      index.re <- Fa>=1|Fa<=0
      geno <- geno[,!index.re]

      #Calculate G matrix
      M <- geno
      p <- ncol(geno)

      cmean<-colMeans(M, na.rm = TRUE)
      V <- apply(M,2,var)
      W <- sweep(sweep(M,2,cmean),2,V, "/")
      W[is.na(W)] <- 0
      W <- as.matrix(W)

      GEMMA2 <- tcrossprod(W)/p

      return(GEMMA2)}

    GAB_matrix<-function(geno,maf){
      #--------------------------------------------------------------------------------------------------------
      # Object: To calculate the GAB matrix
      # Method source:
      #	1.Astle W, Balding DJ. Population Structure and Cryptic Relatedness in Genetic Association Studies, Statistical Science 2009;24:451-471.
      # 2.Leutenegger AL, Prum, B. , Emmanuelle Génin, Verny, C. , Lemainque, A. , & Clerget-Darpoux, F. , et al. Estimation of the inbreeding coefficient through use of genomic data, American Journal of Human Genetics 2003;73:0-523.
      #
      # Input:
      # geno: genotype in numeric format, pure 0, 1, 2 matrix; geno is n individual rows by m snps columns
      # maf: a vector; maf is allele frequency of second allele
      #
      # Output:
      # GAB matrix
      #--------------------------------------------------------------------------------------------------------
      #Remove invariants
      Fa <- colMeans(geno,na.rm=T)/2
      index.re <- Fa>=1|Fa<=0
      geno <- geno[,!index.re]
      maf <- colMeans(geno,na.rm=T)/2

      #Calculate G matrix
      W <- geno
      nSNP <- ncol(W)

      Adj <- .5*maf*(2-maf)
      W <- sweep(sweep(W,2,maf),2,sqrt(Adj), "/")
      W[is.na(W)] <- 0
      W <- as.matrix(W)

      GAB <- tcrossprod(W)/nSNP

      return(GAB)}

    GeneticDistance_matrix<-function(geno){
      #--------------------------------------------------------------------------------------------------------
      # Object: To calculate the Genetic Distance matrix
      # Method source:
      #	Nei M. Genetic Distance between Populations, American Naturalist 1972;106:283-292.
      #
      # Input:
      # geno: genotype in numeric format, pure 0, 1, 2 matrix; geno is n individual rows by m snps columns
      #
      # Output:
      # Genetic Distance matrix
      #--------------------------------------------------------------------------------------------------------
      #Remove invariants
      Fa <- colMeans(geno,na.rm=T)/2
      index.re <- Fa>=1|Fa<=0
      geno <- geno[,!index.re]

      #Calculate G matrix
      marker <- scale(geno,center=TRUE,scale=TRUE)
      Dist <- (as.matrix(dist(marker, method='euclidean'))**2)/ncol(marker)

      GeneticDistance <- exp(-1*Dist)

      return(GeneticDistance)}

    GIBS_matrix<-function(geno){
      #--------------------------------------------------------------------------------------------------------
      # Object: To calculate the GIBS matrix
      #
      # Input:
      # geno: genotype in numeric format, pure 0, 1, 2 matrix; geno is n individual rows by m snps columns
      #
      # Output:
      # GIBS matrix
      #--------------------------------------------------------------------------------------------------------
      #Remove invariants
      Fa <- colMeans(geno,na.rm=T)/2
      index.re <- Fa>=1|Fa<=0
      geno <- geno[,!index.re]

      #Define IBS matrix function
      IBS.matrix <- function(x){
        N = nrow(x)
        Xt<- matrix(NA,nrow= N, ncol=N)
        if (is.matrix(x)!=TRUE) {XX <- as.matrix(x)}
        else {XX<-x}
        for (i in 1:(N-1)) {
          for (j in (i+1):N) {
            Xt[i,j] <- sum(2-abs(XX[i, ] - XX[j, ]))/(2*ncol(XX))
          }
        }
        tt <- t(Xt)
        tt[upper.tri(tt)] <- Xt[upper.tri(Xt)]
        diag(tt)<-c(rep(1,N))
        colnames(tt)<-rownames(x)
        rownames(tt)<-rownames(x)
        return(tt)
      }

      #Calculate G matrix
      IBS <- IBS.matrix(x=geno)
      GIBS <- 2*IBS

      return(GIBS)}

    KING_matrix<-function(geno,maf){
      #--------------------------------------------------------------------------------------------------------
      # Object: To calculate the KING matrix
      # Method source:
      #	Manichaikul A, Mychaleckyj JC, Rich SS et al. Robust relationship inference in genome-wide association studies, Bioinformatics 2010;26:2867-2873.
      #
      # Input:
      # geno: genotype in numeric format, pure 0, 1, 2 matrix; geno is n individual rows by m snps columns
      # maf: a vector; maf is allele frequency of second allele
      #
      # Output:
      # KING matrix
      #--------------------------------------------------------------------------------------------------------
      #Remove invariants
      Fa <- colMeans(geno,na.rm=T)/2
      index.re <- Fa>=1|Fa<=0
      geno <- geno[,!index.re]
      maf <- colMeans(geno,na.rm=T)/2

      #Calculate G matrix
      Adj <- 4*sum(maf*(1-maf))
      N <- nrow(geno)
      Xt <- matrix(NA,nrow= N, ncol=N)

      for (i in 1:(N-1)) {
        for (j in (i+1):N) {
          Xt[i,j] <- 0.5-sum((geno[i, ] - geno[j, ])^2)/Adj
        }
      }

      tt <- t(Xt)
      tt[upper.tri(tt)] <- Xt[upper.tri(Xt)]
      diag(tt)<-c(rep(1,N))
      colnames(tt)<-rownames(geno)
      rownames(tt)<-rownames(geno)

      KING<-tt

      return(KING)}

    HIERFSTAT_matrix<-function(geno){
      #--------------------------------------------------------------------------------------------------------
      # Object: To calculate the HIERFSTAT matrix
      # Method source:
      #	1.Weir BS, Goudet J. A Unified Characterization of Population Structure and Relatedness, Genetics 2017;206:2085-2103.
      #   2.Goudet J. hierfstat, a package for R to compute and test hierarchical F‐statistics, Molecular Ecology Notes 2005;5:184-186.
      #
      # Input:
      #   geno: genotype in numeric format, pure 0, 1, 2 matrix; geno is n individual rows by m snps columns
      #
      # Output:
      # HIERFSTAT matrix
      #--------------------------------------------------------------------------------------------------------
      #Remove invariants
      Fa <- colMeans(geno,na.rm=T)/2
      index.re <- Fa>=1|Fa<=0
      geno <- geno[,!index.re]

      dosage <- geno
      nl <- ncol(dosage)

      #uses matching proba -same equation as for population i.e. Mij=[xiXj+(2-xi)(2-xj)]/4
      Mij <- (tcrossprod(dosage)+tcrossprod(2-dosage))/4
      Ms <- mean(Mij,na.rm=T)

      #Calculate G matrix
      HIERFSTAT <- (Mij-Ms)/(nl-Ms)

      return(HIERFSTAT)}

    haplotype_matrix<-function(geno){
      #--------------------------------------------------------------------------------------------------------
      # Object: To calculate the one-SNP haplotype matrix
      # Method source:
      #	Ferdosi MH, Henshall J, Tier B. Study of the optimum haplotype length to build genomic relationship matrices, Genet Sel Evol 2016;48:75.
      #
      # Input:
      # geno: genotype in numeric format, pure 0, 1, 2 matrix; geno is n individual rows by m snps columns
      #
      # Output:
      # one-SNP haplotype matrix
      #--------------------------------------------------------------------------------------------------------

      #Remove invariants
      Fa <- colMeans(geno,na.rm=T)/2
      index.re <- Fa>=1|Fa<=0
      geno <- geno[,!index.re]

      M <- geno
      nInd <- nrow(M)
      nSNP <- ncol(M)

      #Calculate haplotype matrix
      Haplotype <- hsphase::aio(M)
      XH <- Haplotype

      #Calculate K matrix
      KR <- matrix(rep(1,2),1,2)
      Iah <- matlab::eye(nInd,nInd)
      K <- kronecker(Iah,KR)

      #Calculate T matrix
      Jhm <- matrix(1,2*nInd,nSNP)
      T <- (tcrossprod(Jhm,Jhm)+XH%*%t(XH-Jhm)+t(XH%*%t(XH-Jhm)))/nSNP

      #Calculate G matrix
      L <- K%*%T
      haplotype <- tcrossprod(L,K)

      return(haplotype)}

    GAPIT1_matrix<-function(geno){
      #--------------------------------------------------------------------------------------------------------
      # Object: To calculate the GAPIT1 matrix
      # Method source:
      #	Tang Y, Liu X, Wang J et al. GAPIT Version 2: An Enhanced Integrated Tool for Genomic Association and Prediction, Plant Genome 2016;9.
      #
      # Input:
      # geno: genotype in numeric format, pure 0, 1, 2 matrix; geno is n individual rows by m snps columns
      #
      # Output:
      # GAPIT1 matrix
      #--------------------------------------------------------------------------------------------------------
      #Remove invariants
      Fa <- colMeans(geno,na.rm=T)/2
      index.re <- Fa>=1|Fa<=0
      geno <- geno[,!index.re]

      #Calculate inbreeding coefficient
      het <- 1-abs(geno-1)
      ind.sum <- rowSums(het)
      fi <- ind.sum/(2*ncol(geno))
      inbreeding <- 1-min(fi)

      #Calculate initial G matrix
      nSNP <- ncol(geno)
      nInd <- nrow(geno)
      n <- nInd
      snpMean <- apply(geno,2,mean)   #get mean for each snp
      geno <- t(geno)-snpMean    #operation on matrix and vector goes in direction of column
      G <- crossprod((geno), (geno))

      #Extract diagonals
      i <- 1:n
      j <- (i-1)*n
      index <- i+j
      d <- G[index]
      DL <- min(d)
      DU <- max(d)
      floor <- min(G)

      #Set range between 0 and 2
      top <- 1+inbreeding
      G <- top*(G-floor)/(DU-floor)
      Dmin <- top*(DL-floor)/(DU-floor)

      #Adjust based on expected minimum diagonal
      if(Dmin<1) {
        G[index] <- (G[index]-Dmin+1)/((top+1-Dmin)*.5)
        G[-index] <- G[-index]*(1/Dmin)
      }

      #Limiting the maximum offdiagonal to the top
      Omax <- max(G[-index])
      if(Omax>top) {
        G[-index]=G[-index]*(top/Omax)
      }

      GAPIT1<-G

      return(GAPIT1)}

    GAPIT2_matrix<-function(geno,EV=NULL,nPCs=0){
      #--------------------------------------------------------------------------------------------------------
      # Object: To calculate the GAPIT2 matrix
      # Method source:
      #	Tang Y, Liu X, Wang J et al. GAPIT Version 2: An Enhanced Integrated Tool for Genomic Association and Prediction, Plant Genome 2016;9.
      #
      # Input:
      # geno: genotype in numeric format, pure 0, 1, 2 matrix; geno is n individual rows by m snps columns
      # nPCs: the number of front PCs excluded to calculate kinship
      # EV: Eigen values
      #
      # Output:
      # GAPIT2 matrix
      #--------------------------------------------------------------------------------------------------------
      #get PCA and EV
      PCA <- prcomp(geno,scale.=TRUE)
      PCs <- PCA$x
      Total.number.PCs <- ncol(PCs)
      n <- nrow(PCs)

      #Choose Total.number.PCs-nPCs PCs and EV to calculate G
      sep.PCs <- PCs[,(nPCs+1):(Total.number.PCs)]
      sep.EV <- EV[(nPCs+1):Total.number.PCs]
      Weighted.sep.EV <- sep.EV/sum(sep.EV)

      #X=t(t(sep.PCs)*Weighted.sep.EV)
      X <- sep.PCs
      XMean <- apply(X,2,mean)
      X <- as.matrix(X-XMean)
      G <- tcrossprod((X), (X))

      #Extract diagonals
      i <- 1:n
      j <- (i-1)*n
      index <- i+j
      d <- G[index]
      DL <- min(d)
      DU <- max(d)
      floor <- min(G)
      G <- (G-floor)/(DL-floor)
      MD <- (DU-floor)/(DL-floor)
      if(MD>2) {
        G[index]=G[index]/(MD-1)+1
      }

      GAPIT2<-G

      return (GAPIT2)}


    if (method=="VanRaden1"){
      G<-VanRaden1_matrix(geno,maf)}
    if (method=="VanRaden2"){
      G<-VanRaden2_matrix(geno)}
    if (method=="VanRaden3"){
      G<-VanRaden3_matrix(geno,A)}
    if (method=="GEMMA1"){
      G<-GEMMA1_matrix(geno)}
    if (method=="GEMMA2"){
      G<-GEMMA2_matrix(geno)}
    if (method=="GAB"){
      G<-GAB_matrix(geno,maf)}
    if (method=="Genetic Distance"){
      G<-GeneticDistance_matrix(geno)}
    if (method=="GIBS"){
      G<-GIBS_matrix(geno)}
    if (method=="KING"){
      G<-KING_matrix(geno,maf)}
    if (method=="Weighted"){
      G<-Weighted_matrix(geno,maf)}
    if (method=="HIERFSTAT"){
      G<-HIERFSTAT_matrix(geno)}
    if (method=="one-SNP haplotype"){
      G<-haplotype_matrix(geno)}
    if (method=="GAPIT1"){
      G<-GAPIT1_matrix(geno)}
    if (method=="GAPIT2"){
      G<-GAPIT2_matrix(geno)}


    return(G)}
