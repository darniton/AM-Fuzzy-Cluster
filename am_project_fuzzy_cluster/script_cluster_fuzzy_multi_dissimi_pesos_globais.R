#Importing data into the project

fac <- read.csv(file="dataset/mfeat_fac.txt", header = T, sep = ";")
fou <- read.csv(file="dataset/mfeat_fou.txt", header = T, sep = ";")
kar <- read.csv(file="dataset/mfeat_kar.txt", header = T, sep = ";")

#Calculating the dissimilarity matrix
m_fac <- as.matrix( dist(fac[,1:216], method = "euclidean"))
m_fou <- as.matrix( dist(fou[,1:76], method = "euclidean"))
m_kar <- as.matrix( dist(kar[,1:64], method = "euclidean"))


#Initialization

#Number of clusters
K <- 10 

#parameter that controls the fuzziness of membership for each object
m <- 1.6

#Interaction threshold
t <- 300

#prototype cardinality
q <- 3

#Number of objects
numeroObjetos <- 2000

#value of s
s <- 1

#set lambada
#lambda <- c(1,1,1)

#Selecting randomly K different prototypes
gerarPrototiposIniciais <- function( numeroCluster = K, card = q, objetos = numeroObjetos) {
  
  matrizG = matrix(sample(1:objetos, numeroCluster*card, replace=FALSE), ncol = card, byrow = T)
  matrizG
}

#For each object, compute the degree of pertinence for each cluster Ck
gerarMatrizU <- function( numeroCluster = K, objetos = numeroObjetos, matrizPrototipo , matrizDissimilaridade01, matrizDissimilaridade02, matrizDissimilaridade03, matrizLambida ) {
  
  #initializes the matrix U
  Ui <- c()
  
  #calculates uik for each cluster
  for ( cluster in seq(1,numeroCluster)) {
    
    #computes the ui of a cluster
    for (objs in seq(1,objetos)) {
      
      #calculating the distance of each object in relation to the representatives of the cluster
      u <- ( (matrizLambida[1] * sum(matrizDissimilaridade01[objs , matrizPrototipo [cluster,1]], matrizDissimilaridade01[objs , matrizPrototipo [cluster,2]], matrizDissimilaridade01[objs , matrizPrototipo[cluster,3]])) +
               (matrizLambida[2] * sum(matrizDissimilaridade02[objs , matrizPrototipo [cluster,1]], matrizDissimilaridade02[objs , matrizPrototipo [cluster,2]], matrizDissimilaridade02[objs , matrizPrototipo[cluster,3]])) +
               (matrizLambida[3] * sum(matrizDissimilaridade03[objs , matrizPrototipo [cluster,1]], matrizDissimilaridade03[objs , matrizPrototipo [cluster,2]], matrizDissimilaridade03[objs , matrizPrototipo[cluster,3]])) )
      
      total <- 0
      ut <- 0
      
      #computes the distance of for each h ranging from 1 to k
      for (clusterTotal in seq(1,numeroCluster)) {
        
        #calculating the distance of each object in relation to each cluster individually
        total <- ( (matrizLambida[1] * sum(matrizDissimilaridade01[objs,matrizPrototipo[clusterTotal,1]], matrizDissimilaridade01[objs,matrizPrototipo[clusterTotal,2]], matrizDissimilaridade01[objs,matrizPrototipo[clusterTotal,3]])) +
                     (matrizLambida[2] * sum(matrizDissimilaridade02[objs,matrizPrototipo[clusterTotal,1]], matrizDissimilaridade02[objs,matrizPrototipo[clusterTotal,2]], matrizDissimilaridade02[objs,matrizPrototipo[clusterTotal,3]])) +
                     (matrizLambida[3] * sum(matrizDissimilaridade03[objs,matrizPrototipo[clusterTotal,1]], matrizDissimilaridade03[objs,matrizPrototipo[clusterTotal,2]], matrizDissimilaridade03[objs,matrizPrototipo[clusterTotal,3]])) )
        
        #calculation that represents the division of the distance of an object to a specific cluster by the distance of this object to each cluster individually
        ut <- sum(ut, (u/total)^(1/(m-1)))
        
      }
      
      #calculates the inverse for the purpose of specifying that the longest distance represents the least degree of pertinence
      Ui <- c(Ui, ut^(-1))
    }
    
    
  }
  
  matrizUi <- matrix(Ui, ncol = numeroCluster)
  return(matrizUi)
}

#Objective function

gerarFuncaoObjetivo <- function( numeroCluster = K, objetos = numeroObjetos,  matrizDissimilaridade01, matrizDissimilaridade02, matrizDissimilaridade03, matrizLambida, matrizPrototipo, matrizu) {
  
  j <- 0
  
  
  for (cluster in seq(1,numeroCluster)) {
    
    gu <- 0
    
    
    for (objs in seq(1,objetos)) {
      
      #performs the sum of the product of each element of the membership matrix by the distance of each object to a cluster
      gu <- sum(gu, ((matrizu[objs, cluster])^m * matrizLambida[1] * (sum(matrizDissimilaridade01[objs,matrizPrototipo[cluster,1]], matrizDissimilaridade01[objs,matrizPrototipo[cluster,2]], matrizDissimilaridade01[objs,matrizPrototipo[cluster,3]]))) ,
                ((matrizu[objs, cluster])^m * matrizLambida[2] * (sum(matrizDissimilaridade02[objs,matrizPrototipo[cluster,1]], matrizDissimilaridade02[objs,matrizPrototipo[cluster,2]], matrizDissimilaridade02[objs,matrizPrototipo[cluster,3]]))) ,
                ((matrizu[objs, cluster])^m * matrizLambida[3] * (sum(matrizDissimilaridade03[objs,matrizPrototipo[cluster,1]], matrizDissimilaridade03[objs,matrizPrototipo[cluster,2]], matrizDissimilaridade03[objs,matrizPrototipo[cluster,3]]))) )
      
      
      
    }
    
    j <- sum(j, gu)
    
  }
  
  return(j)
}

#Otimizacao

#Computing the best prototype

computarMelhorPrototipo <- function(matrizu, matrizDissimilaridade01, matrizDissimilaridade02, matrizDissimilaridade03, cluster, objetos, matrizLambida){
  
  s <- c()
  
  for (objs in seq(1:objetos)) {
    
    s <- c(s, sum( (((matrizu[,cluster])^m) * matrizLambida[1] * matrizDissimilaridade01[1:objetos,objs]),
                   (((matrizu[,cluster])^m) * matrizLambida[2] * matrizDissimilaridade02[1:objetos,objs]), 
                   (((matrizu[,cluster])^m) * matrizLambida[3] * matrizDissimilaridade03[1:objetos,objs])  ))
    
    
  }
  
  return(order(s)[1:3])
  
  
  
}


#Better prototypes

gerarPrototiposMelhorados <- function(numeroCluster = K, objetos = numeroObjetos,  matrizDissimilaridade01, matrizDissimilaridade02, matrizDissimilaridade03, matrizu, matrizLambida, card = q) {
  
  matrizGMelhorada <- c()
  matrizGM <- c()
  for (clust in seq(1,numeroCluster)){
    
    matrizGM <- c(matrizGM, computarMelhorPrototipo(matrizDissimilaridade01 = matrizDissimilaridade01, matrizDissimilaridade02 = matrizDissimilaridade02, matrizDissimilaridade03 =matrizDissimilaridade03, matrizLambida = matrizLambida, matrizu = matrizu, objetos = objetos, cluster = clust))
    
  }
  
  matrizGMelhorada <- matrix(matrizGM, ncol = card, byrow = T)
  return(matrizGMelhorada)
  
}

#calculating the weight matrix

gerarMatrizLambda <- function( numeroCluster = K, objetos = numeroObjetos, matrizPrototipo , matrizDissimilaridade01, matrizDissimilaridade02, matrizDissimilaridade03, matrizu ) {
  
  #initializes the matrix U
  li <- c()
  p1 <- 0
  p2 <- 0
  p3 <- 0
  
  #computes the lkj for each cluster
  for ( cluster in seq(1,numeroCluster)) {
    
    #computes the li of a cluster
    for (objs in seq(1,objetos)) {
      
      #calculating the distance of each object in relation to the representatives of the cluster
      p1 <- sum(p1, ((matrizu[objs, cluster]^m) * sum(matrizDissimilaridade01[objs , matrizPrototipo [cluster,1]], matrizDissimilaridade01[objs , matrizPrototipo [cluster,2]], matrizDissimilaridade01[objs , matrizPrototipo[cluster,3]])))
      p2 <- sum(p2, ((matrizu[objs, cluster]^m) * sum(matrizDissimilaridade02[objs , matrizPrototipo [cluster,1]], matrizDissimilaridade02[objs , matrizPrototipo [cluster,2]], matrizDissimilaridade02[objs , matrizPrototipo[cluster,3]])))
      p3 <- sum(p3, ((matrizu[objs, cluster]^m) * sum(matrizDissimilaridade03[objs , matrizPrototipo [cluster,1]], matrizDissimilaridade03[objs , matrizPrototipo [cluster,2]], matrizDissimilaridade03[objs , matrizPrototipo[cluster,3]])))
      
    }  
    
    
  }
  
  p <- (p1*p2*p3)^(1/3)
  matrizLi <- c(p/p1, p/p2, p/p3)
  return(matrizLi)
}



gerarCluster <- function( nint = t, numeroCluster = K, objetos = numeroObjetos,  matrizDissimilaridade01 = m_fac, matrizDissimilaridade02 = m_fou, matrizDissimilaridade03 = m_kar, e = 0.01 ){
  
  
  #inicializacao
  
  L0 <- c(1,1,1)
  G0 <- gerarPrototiposIniciais(numeroCluster = numeroCluster, card = 3, objetos = objetos)
  U0 <- gerarMatrizU(numeroCluster = numeroCluster, objetos = objetos, matrizPrototipo = G0 , matrizDissimilaridade01 = matrizDissimilaridade01, matrizDissimilaridade02 = matrizDissimilaridade02, matrizDissimilaridade03 = matrizDissimilaridade03, matrizLambida = L0)
  J0 <- gerarFuncaoObjetivo(numeroCluster = numeroCluster, objetos = objetos,  matrizDissimilaridade01 = matrizDissimilaridade01, matrizDissimilaridade02 = matrizDissimilaridade02, matrizDissimilaridade03 = matrizDissimilaridade03, matrizLambida = L0, matrizPrototipo = G0, matrizu = U0)
  
  L <- L0
  G <- G0
  U <- U0
  J <- J0
  
  np <- 0
  Jt <- 0
  Jt1 <- 0
  
  repeat {
    
    np <- sum(np,1)
    
    if (np == 1 ){
      
      
      Gt <- gerarPrototiposMelhorados(numeroCluster = numeroCluster, objetos = objetos,  matrizDissimilaridade01 = matrizDissimilaridade01, matrizDissimilaridade02 = matrizDissimilaridade02, matrizDissimilaridade03 = matrizDissimilaridade03, matrizu = U0, matrizLambida = L0, card = 3)
      Lt <- gerarMatrizLambda(numeroCluster = numeroCluster, objetos = objetos, matrizPrototipo = Gt , matrizDissimilaridade01 = matrizDissimilaridade01, matrizDissimilaridade02 = matrizDissimilaridade02, matrizDissimilaridade03 = matrizDissimilaridade03, matrizu = U0 )
      Ut <- gerarMatrizU(numeroCluster = numeroCluster, objetos = objetos, matrizPrototipo = Gt , matrizDissimilaridade01 = matrizDissimilaridade01, matrizDissimilaridade02 = matrizDissimilaridade02, matrizDissimilaridade03 = matrizDissimilaridade03, matrizLambida = Lt)
      Jt <- gerarFuncaoObjetivo(numeroCluster = numeroCluster, objetos = objetos,  matrizDissimilaridade01 = matrizDissimilaridade01, matrizDissimilaridade02 = matrizDissimilaridade02, matrizDissimilaridade03 = matrizDissimilaridade03, matrizPrototipo = Gt, matrizu = Ut, matrizLambida = Lt)
      
      L <- Lt
      G <- Gt
      U <- Ut
      J <- Jt
      
    } else {
      
      
      if (np%%2 == 0) {
        
        Gt1 <- gerarPrototiposMelhorados(numeroCluster = numeroCluster, objetos = objetos,  matrizDissimilaridade01 = matrizDissimilaridade01, matrizDissimilaridade02 = matrizDissimilaridade02, matrizDissimilaridade03 = matrizDissimilaridade03, matrizu = Ut, matrizLambida = Lt , card = 3)
        Lt1 <- gerarMatrizLambda(numeroCluster = numeroCluster, objetos = objetos, matrizPrototipo = Gt1 , matrizDissimilaridade01 = matrizDissimilaridade01, matrizDissimilaridade02 = matrizDissimilaridade02, matrizDissimilaridade03 = matrizDissimilaridade03, matrizu = Ut )
        Ut1 <- gerarMatrizU(numeroCluster = numeroCluster, objetos = objetos, matrizPrototipo = Gt1 , matrizDissimilaridade01 = matrizDissimilaridade01, matrizDissimilaridade02 = matrizDissimilaridade02, matrizDissimilaridade03 = matrizDissimilaridade03, matrizLambida = Lt1)
        Jt1 <- gerarFuncaoObjetivo(numeroCluster = numeroCluster, objetos = objetos,  matrizDissimilaridade01 = matrizDissimilaridade01, matrizDissimilaridade02 = matrizDissimilaridade02, matrizDissimilaridade03 = matrizDissimilaridade03, matrizPrototipo = Gt1, matrizu = Ut1, matrizLambida = Lt1)
        
        L <- Lt1
        G <- Gt1
        U <- Ut1
        J <- Jt1
        
      } else {
        
        
        Gt <- gerarPrototiposMelhorados(numeroCluster = numeroCluster, objetos = objetos,  matrizDissimilaridade01 = matrizDissimilaridade01, matrizDissimilaridade02 = matrizDissimilaridade02, matrizDissimilaridade03 = matrizDissimilaridade03, matrizu = Ut1, matrizLambida = Lt1, card = 3)
        Lt <- gerarMatrizLambda(numeroCluster = numeroCluster, objetos = objetos, matrizPrototipo = Gt , matrizDissimilaridade01 = matrizDissimilaridade01, matrizDissimilaridade02 = matrizDissimilaridade02, matrizDissimilaridade03 = matrizDissimilaridade03, matrizu = Ut1 )
        Ut <- gerarMatrizU(numeroCluster = numeroCluster, objetos = objetos, matrizPrototipo = Gt , matrizDissimilaridade01 = matrizDissimilaridade01, matrizDissimilaridade02 = matrizDissimilaridade02, matrizDissimilaridade03 = matrizDissimilaridade03, matrizLambida = Lt)
        Jt <- gerarFuncaoObjetivo(numeroCluster = numeroCluster, objetos = objetos,  matrizDissimilaridade01 = matrizDissimilaridade01, matrizDissimilaridade02 = matrizDissimilaridade02, matrizDissimilaridade03 = matrizDissimilaridade03, matrizPrototipo = Gt, matrizu = Ut, matrizLambida = Lt)
        
        L <- Lt
        G <- Gt
        U <- Ut
        J <- Jt
        
      }
      
      
    }
    
    
    if (abs(Jt1-Jt) <= e | np == nint) break()
    
    
  }
  
  return(result <- (list(L, G, U, J)))
  
}


MFCMdd_RWG_P <- function(nrep, nint = t, numeroCluster = K, objetos = numeroObjetos,  matrizDissimilaridade01 = m_fac, matrizDissimilaridade02 = m_fou, matrizDissimilaridade03 = m_kar, e = 0.01) {
  
  
  for (i in seq(1:nrep)) { 
    
    if (i == 1) {
    
      resultadoParcial <- gerarCluster(nint = nint, numeroCluster = numeroCluster, objetos = objetos,  matrizDissimilaridade01 = matrizDissimilaridade01, matrizDissimilaridade02 = matrizDissimilaridade02, matrizDissimilaridade03 = matrizDissimilaridade03, e = e) 
      result <- resultadoParcial
      
    } else{
      
      resultadoParcial <- gerarCluster(nint = nint, numeroCluster = numeroCluster, objetos = objetos,  matrizDissimilaridade01 = matrizDissimilaridade01, matrizDissimilaridade02 = matrizDissimilaridade02, matrizDissimilaridade03 = matrizDissimilaridade03, e = e) 
      result <- c(result, resultadoParcial)
      
    }
  }

  resultFuncaObjetivo <- c()
  
  
  for (j in seq(1:nrep)){
    
    
    resultFuncaObjetivo <- c(resultFuncaObjetivo, result[[j*4]])
    
    
  }
  
  contador <- c(0:nrep)
  
  if (order(resultFuncaObjetivo)[1] == 1 ) {
    
    L <- result[[ 1 ]] 
    G <- result[[ 2 ]]
    U <- result[[ 3 ]] 
    J <- result[[ 4 ]]
    
  } else {
    
    L <- result[[ (((order(resultFuncaObjetivo)[1])-1)*4)+1 ]]
    G <- result[[ (((order(resultFuncaObjetivo)[1])-1)*4)+2 ]]
    U <- result[[ (((order(resultFuncaObjetivo)[1])-1)*4)+3 ]]
    J <- result[[ (((order(resultFuncaObjetivo)[1])-1)*4)+4 ]]
    
  }
  
  clusterHard <- c()
   
  for (i in seq(1:objetos)){
    
    clusterHard <- c(clusterHard, c(i, order(U[i,])[numeroCluster]))
    
  }
  
  clusterHardMatriz <- matrix(clusterHard, ncol = 2, nrow =objetos, byrow = T )
  
  return(list(L, G, U, J, clusterHardMatriz))
  
}


