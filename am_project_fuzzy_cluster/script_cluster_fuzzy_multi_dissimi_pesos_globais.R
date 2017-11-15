library(clues)

#Initialization

#Number of clusters
k <- 10 

#parameter that controls the fuzziness of membership for each object
M <- 1.6

#Interaction threshold
T <- 100

#prototype cardinality
q <- 3

#Number of objects
numeroObjetos <- 2000

#value of s
s <- 1

#set lambada
lambda <- c(1,1,1)

#epsilon
E = 10^(-100)

#interaction
t <- 20


#Importing data into the project
fac <- read.csv(file="mfeat_fac_tratado.txt", header = T, sep = ";")
fou <- read.csv(file="mfeat_fou_tratado.txt", header = T, sep = ";")
kar <- read.csv(file="mfeat_kar_tratado.txt", header = T, sep = ";")

#class label
rotulo <- fac[,218]

#data without the class label
fac <- fac[,1:216]
fou <- fou[,1:76]
kar <- kar[,1:64]

#normalization
range01 <- function(x){(x-min(x))/(max(x)-min(x))}

fac <- range01(fac)
fou <- range01(fou)
kar <- range01(kar)

#
normMatrizDissi <- function(matrizDissimilaridade, objetos = 2000  ) {
  
  s <- c()
  
  for (objs in seq(1:objetos)) {
    
    s <- c(s, sum( matrizDissimilaridade[,objs])  )
    
  }
  
  #calculo de T
  
  T <- sum(matrizDissimilaridade[order(s)[1],])
  
  return(T)
  
}

m_fac <- as.matrix( dist(fac[,1:216], method = "euclidean"))
m_fou <- as.matrix( dist(fou[,1:76], method =  "euclidean"))
m_kar <- as.matrix( dist(kar[,1:64], method = "euclidean"))

m_fac <- m_fac / normMatrizDissi(matrizDissimilaridade = m_fac)
m_fou <- m_fou / normMatrizDissi(matrizDissimilaridade = m_fou)
m_kar <- m_kar / normMatrizDissi(matrizDissimilaridade = m_kar)


#Selecting randomly K different prototypes
gerarPrototiposIniciais <- function( numeroCluster = k, card = q, objetos = numeroObjetos) {
  
  matrizG = matrix(sample(1:objetos, numeroCluster*card, replace=FALSE), ncol = card, byrow = T)
  return(matrizG)
}

#For each object, compute the degree of pertinence for each cluster Ck
gerarMatrizU <- function( numeroCluster = k, objetos = numeroObjetos, matrizPrototipo , matrizDissimilaridade01 = m_fac, matrizDissimilaridade02 = m_fou, matrizDissimilaridade03 = m_kar, matrizLambida , m = M) {
  
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

gerarFuncaoObjetivo <- function( numeroCluster = k, objetos = numeroObjetos,  matrizDissimilaridade01 = m_fac, matrizDissimilaridade02 = m_fou, matrizDissimilaridade03 = m_kar, matrizLambida, matrizPrototipo, matrizu, m = M) {
  
  j <- 0
  
  
  for (cluster in seq(1,numeroCluster)) {
    
    gu <- 0
    
    
    for (objs in seq(1,objetos)) {
      
      #performs the sum of the product of each element of the membership matrix by the distance of each object to a cluster
      gu <- sum(gu, sum( ((matrizu[objs, cluster])^m * matrizLambida[1] * (sum(matrizDissimilaridade01[objs,matrizPrototipo[cluster,1]], matrizDissimilaridade01[objs,matrizPrototipo[cluster,2]], matrizDissimilaridade01[objs,matrizPrototipo[cluster,3]]))) ,
                         ((matrizu[objs, cluster])^m * matrizLambida[2] * (sum(matrizDissimilaridade02[objs,matrizPrototipo[cluster,1]], matrizDissimilaridade02[objs,matrizPrototipo[cluster,2]], matrizDissimilaridade02[objs,matrizPrototipo[cluster,3]]))) ,
                         ((matrizu[objs, cluster])^m * matrizLambida[3] * (sum(matrizDissimilaridade03[objs,matrizPrototipo[cluster,1]], matrizDissimilaridade03[objs,matrizPrototipo[cluster,2]], matrizDissimilaridade03[objs,matrizPrototipo[cluster,3]]))) )  )
      
      
      
    }
    
    j <- sum(j, gu)
    
  }
  
  #print(j)
  return(j)
}

#Otimizacao

#Computing the best prototype

#computarMelhorPrototipo <- function(matrizu, matrizDissimilaridade01 = m_fac, matrizDissimilaridade02= m_fou, matrizDissimilaridade03=m_kar, cluster=0 , objetos = 2000, matrizLambida, card = 3, m = 1.6){

#  s <- c()

#print(matrizu[1:objetos,cluster]^m)

#  for (objs in seq(1:objetos)) {



#    s <- c(s, sum( (matrizu[1:objetos,cluster]^m) * ( matrizLambida[1] * matrizDissimilaridade01[1:objetos,objs] +
#                                                      matrizLambida[2] * matrizDissimilaridade02[1:objetos,objs] + 
#                                                      matrizLambida[3] * matrizDissimilaridade03[1:objetos,objs]  )  )  )


#  }

#  return(order(s)[1:card])



#}


computarMelhorPrototipo <- function(matrizu, matrizLambida, cluster, card = q, objetos = numeroObjetos, matrizDissimilaridade01 = m_fac, matrizDissimilaridade02 = m_fou, matrizDissimilaridade03 = m_kar, m = M) {
  
  u <- (matrizu[,cluster])^m
  z <- matrizLambida[1]*matrizDissimilaridade01 + matrizLambida[2]*matrizDissimilaridade02 + matrizLambida[3]*matrizDissimilaridade03
  l <- u*z
  g <- c()
  
  for (i in seq(1:objetos)) {
    
    g <- c(g, sum(l[,i]))
    
  }
  
  return(order(g)[1:card])
  
}



#Better prototypes

gerarPrototiposMelhorados <- function(numeroCluster = k, objetos = numeroObjetos,  matrizDissimilaridade01 = m_fac, matrizDissimilaridade02 = m_fou, matrizDissimilaridade03 = m_kar, matrizu, matrizLambida, card = q, m = M) {
  
  matrizGMelhorada <- c()
  matrizGM <- c()
  for (clust in seq(1,numeroCluster)){
    
    matrizGM <- c(matrizGM, computarMelhorPrototipo(matrizDissimilaridade01 = matrizDissimilaridade01, matrizDissimilaridade02 = matrizDissimilaridade02, matrizDissimilaridade03 =matrizDissimilaridade03, matrizLambida = matrizLambida, matrizu = matrizu, objetos = objetos, cluster = clust , card = card, m = m))
    
  }
  
  matrizGMelhorada <- matrix(matrizGM, ncol = card, byrow = T)
  return(matrizGMelhorada)
  
}

#calculating the weight matrix

gerarMatrizLambda <- function( numeroCluster = k, objetos = numeroObjetos, matrizPrototipo , matrizDissimilaridade01=m_fac, matrizDissimilaridade02=m_fou, matrizDissimilaridade03=m_kar, matrizu, m = M ) {
  
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



gerarCluster <- function( nint = t, numeroCluster = k, objetos = numeroObjetos,  matrizDissimilaridade01 = m_fac, matrizDissimilaridade02 = m_fou, matrizDissimilaridade03 = m_kar, e = E , m = M, card = q ){
  
  
  #inicializacao
  
  L0 <- c(1,1,1)
  G0 <- gerarPrototiposIniciais(numeroCluster = numeroCluster, card = card, objetos = objetos)
  U0 <- gerarMatrizU(numeroCluster = numeroCluster, objetos = objetos, matrizPrototipo = G0 , matrizDissimilaridade01 = matrizDissimilaridade01, matrizDissimilaridade02 = matrizDissimilaridade02, matrizDissimilaridade03 = matrizDissimilaridade03, matrizLambida = L0, m = m)
  J0 <- gerarFuncaoObjetivo(numeroCluster = numeroCluster, objetos = objetos,  matrizDissimilaridade01 = matrizDissimilaridade01, matrizDissimilaridade02 = matrizDissimilaridade02, matrizDissimilaridade03 = matrizDissimilaridade03, matrizLambida = L0, matrizPrototipo = G0, matrizu = U0, m = m)
  
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
      
      
      Gt <- gerarPrototiposMelhorados(numeroCluster = numeroCluster, objetos = objetos,  matrizDissimilaridade01 = matrizDissimilaridade01, matrizDissimilaridade02 = matrizDissimilaridade02, matrizDissimilaridade03 = matrizDissimilaridade03, matrizu = U0, matrizLambida = L0, card = q , m = m )
      Lt <- gerarMatrizLambda(numeroCluster = numeroCluster, objetos = objetos, matrizPrototipo = Gt , matrizDissimilaridade01 = matrizDissimilaridade01, matrizDissimilaridade02 = matrizDissimilaridade02, matrizDissimilaridade03 = matrizDissimilaridade03, matrizu = U0, m = m )
      Ut <- gerarMatrizU(numeroCluster = numeroCluster, objetos = objetos, matrizPrototipo = Gt , matrizDissimilaridade01 = matrizDissimilaridade01, matrizDissimilaridade02 = matrizDissimilaridade02, matrizDissimilaridade03 = matrizDissimilaridade03, matrizLambida = Lt, m = m)
      Jt <- gerarFuncaoObjetivo(numeroCluster = numeroCluster, objetos = objetos,  matrizDissimilaridade01 = matrizDissimilaridade01, matrizDissimilaridade02 = matrizDissimilaridade02, matrizDissimilaridade03 = matrizDissimilaridade03, matrizPrototipo = Gt, matrizu = Ut, matrizLambida = Lt, m = m)
      
      L <- Lt
      G <- Gt
      U <- Ut
      J <- Jt
      
    } else {
      
      
      if (np%%2 == 0) {
        
        Gt1 <- gerarPrototiposMelhorados(numeroCluster = numeroCluster, objetos = objetos,  matrizDissimilaridade01 = matrizDissimilaridade01, matrizDissimilaridade02 = matrizDissimilaridade02, matrizDissimilaridade03 = matrizDissimilaridade03, matrizu = Ut, matrizLambida = Lt , m = M, card = q )
        Lt1 <- gerarMatrizLambda(numeroCluster = numeroCluster, objetos = objetos, matrizPrototipo = Gt1 , matrizDissimilaridade01 = matrizDissimilaridade01, matrizDissimilaridade02 = matrizDissimilaridade02, matrizDissimilaridade03 = matrizDissimilaridade03, matrizu = Ut, m = m )
        Ut1 <- gerarMatrizU(numeroCluster = numeroCluster, objetos = objetos, matrizPrototipo = Gt1 , matrizDissimilaridade01 = matrizDissimilaridade01, matrizDissimilaridade02 = matrizDissimilaridade02, matrizDissimilaridade03 = matrizDissimilaridade03, matrizLambida = Lt1, m = m)
        Jt1 <- gerarFuncaoObjetivo(numeroCluster = numeroCluster, objetos = objetos,  matrizDissimilaridade01 = matrizDissimilaridade01, matrizDissimilaridade02 = matrizDissimilaridade02, matrizDissimilaridade03 = matrizDissimilaridade03, matrizPrototipo = Gt1, matrizu = Ut1, matrizLambida = Lt1 , m = m)
        
        L <- Lt1
        G <- Gt1
        U <- Ut1
        J <- Jt1
        
      } else {
        
        
        Gt <- gerarPrototiposMelhorados(numeroCluster = numeroCluster, objetos = objetos,  matrizDissimilaridade01 = matrizDissimilaridade01, matrizDissimilaridade02 = matrizDissimilaridade02, matrizDissimilaridade03 = matrizDissimilaridade03, matrizu = Ut1, matrizLambida = Lt1, m = M, card = q )
        Lt <- gerarMatrizLambda(numeroCluster = numeroCluster, objetos = objetos, matrizPrototipo = Gt , matrizDissimilaridade01 = matrizDissimilaridade01, matrizDissimilaridade02 = matrizDissimilaridade02, matrizDissimilaridade03 = matrizDissimilaridade03, matrizu = Ut1, m = m )
        Ut <- gerarMatrizU(numeroCluster = numeroCluster, objetos = objetos, matrizPrototipo = Gt , matrizDissimilaridade01 = matrizDissimilaridade01, matrizDissimilaridade02 = matrizDissimilaridade02, matrizDissimilaridade03 = matrizDissimilaridade03, matrizLambida = Lt, m = m)
        Jt <- gerarFuncaoObjetivo(numeroCluster = numeroCluster, objetos = objetos,  matrizDissimilaridade01 = matrizDissimilaridade01, matrizDissimilaridade02 = matrizDissimilaridade02, matrizDissimilaridade03 = matrizDissimilaridade03, matrizPrototipo = Gt, matrizu = Ut, matrizLambida = Lt, m = m)
        
        L <- Lt
        G <- Gt
        U <- Ut
        J <- Jt
        
      }
      
      
    }
    
    
    if (abs(Jt1-Jt) <= e | np == nint) break()
    
    
  }
  
  U <- c()
  return(result <- (list(L, G, U, J)))
  
}


MFCMdd_RWG_P <- function(nrep, nint = t, numeroCluster = k, objetos = numeroObjetos,  matrizDissimilaridade01 = m_fac, matrizDissimilaridade02 = m_fou, matrizDissimilaridade03 = m_kar, e = E, m = M) {
  
  cont <- 0
  for (i in seq(1:nrep)) { 
    
    if (i == 1) {
      
      resultadoParcial <- gerarCluster(nint = nint, numeroCluster = numeroCluster, objetos = objetos,  matrizDissimilaridade01 = matrizDissimilaridade01, matrizDissimilaridade02 = matrizDissimilaridade02, matrizDissimilaridade03 = matrizDissimilaridade03, e = e, m = m) 
      result <- resultadoParcial
      print(resultadoParcial[4])
      print(resultadoParcial[2])
      
      
      
    } else{
      
      resultadoParcial <- gerarCluster(nint = nint, numeroCluster = numeroCluster, objetos = objetos,  matrizDissimilaridade01 = matrizDissimilaridade01, matrizDissimilaridade02 = matrizDissimilaridade02, matrizDissimilaridade03 = matrizDissimilaridade03, e = e, m = m) 
      result <- c(result, resultadoParcial)
      print(resultadoParcial[4])
      print(resultadoParcial[2])
      
    }
    
    cont <- sum(cont, 1)
    print(cont)
    
    
  }
  
  resultFuncaObjetivo <- c()
  
  
  for (j in seq(1:nrep)){
    
    
    resultFuncaObjetivo <- c(resultFuncaObjetivo, result[[j*4]])
    
    
  }
  
  contador <- c(0:nrep)
  
  if (order(resultFuncaObjetivo)[1] == 1 ) {
    
    L <- result[[ 1 ]] 
    G <- result[[ 2 ]]
    U <- gerarMatrizU(numeroCluster = k, objetos = numeroObjetos, matrizPrototipo = G, matrizDissimilaridade01 = m_fac, matrizDissimilaridade02 = m_fou, matrizDissimilaridade03 = m_kar, matrizLambida = L, m = m)
    J <- result[[ 4 ]]
    
  } else {
    
    L <- result[[ (((order(resultFuncaObjetivo)[1])-1)*4)+1 ]]
    G <- result[[ (((order(resultFuncaObjetivo)[1])-1)*4)+2 ]]
    U <- gerarMatrizU(numeroCluster = k, objetos = numeroObjetos, matrizPrototipo = G, matrizDissimilaridade01 = m_fac, matrizDissimilaridade02 = m_fou, matrizDissimilaridade03 = m_kar, matrizLambida = L , m = m)
    J <- result[[ (((order(resultFuncaObjetivo)[1])-1)*4)+4 ]]
    
  }
  
  clusterHard <- c()
  
  for (i in seq(1:objetos)){
    
    clusterHard <- c(clusterHard, c(i, order(U[i,])[numeroCluster]-1))
    
  }
  
  clusterHardMatriz <- matrix(clusterHard, ncol = 2, nrow =objetos, byrow = T )
  rand <- adjustedRand(rotulo, clusterHardMatriz[,2])
  
  return(list(L, G, U, J, clusterHardMatriz, rand))
  #return(list(L, G, U, J))
  
}
