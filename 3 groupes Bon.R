
##### Test Generique de la Fonction #####
MeData <- list()
MeResF <- list()
G <- 3
P <- 6
r <- 2
MatBeta <- matrix(nrow=P,ncol=1)
Rvect <- rep(r/sqrt(P),P)
Beta <- runif(P,-100,100)
MatBeta[,1] <- Beta

for(k in 2:G)
{
  r <- 1
  while(r!=0)
  {
    a <- runif(1,0.5,0.6)
    u <- rbinom(P,1,a)
    u[which(u==0)] <- -1
    Vect_sign_dev <- u
    Dev <- Vect_sign_dev*Rvect
    MatDev <- matrix(rep(Beta+Dev,k-1),ncol=k-1)
    MatTest <- MatBeta-MatDev
    uTest <- apply(MatTest,2,ft)
    r <- length(which(uTest==0))
    if(r==0)
    {
      MatBeta <- cbind(MatBeta,MatDev[,1])
    }
    
  }
}

Vi <- runif(G,0.2,2)
Ve <- runif(G,0.2,2)
Vr <- runif(G,0.2,2)

k <- 99
while(k!=100)
{
  ech <- sample(1:100,500,replace=T)  
  u_ech <- as.vector(table(ech))
  k <- length(u_ech)
}
n_ech <- rep(5,100)+u_ech

Lw <- list(n_ech[1:33],n_ech[34:66],n_ech[67:100])
To <- 10
res <- SiMMMEI(MatBeta,Vi,Vr,Ve,Lw,To)

Data <- res$Data
MeData[[1]] <- res
formula <- YF~V2+V3+V4+V5+V6
ResF <- EstimF(formula,Data,G)
MeResF[[1]] <- ResF


uInit <- c(rep(1,33),rep(2,33),rep(3,34))

Ind_R <- Ind_Regr(ResF$Rg,uInit)
I1 <- Ind_R[1]
I2 <- Ind_R[2]
It <- ResF$Rs
D_Beta <- RMS_Beta(MatBeta[,unique(ResF$Rg)],ResF$theta$MatBeta)
D_Vi <- RMS_Sigma(Vi[unique(ResF$Rg)],ResF$theta$Vi)
D_Vr <- RMS_Sigma(Vr[unique(ResF$Rg)],ResF$theta$Vr)
D_Ve <- RMS_Sigma(Ve[unique(ResF$Rg)],ResF$theta$Ve)

for(kk in 2:100)
{
  G <- 3
  P <- 6
  r <- 2
  MatBeta <- matrix(nrow=P,ncol=1)
  Rvect <- rep(r/sqrt(P),P)
  Beta <- runif(P,-100,100)
  MatBeta[,1] <- Beta
  
  for(k in 2:G)
  {
    r <- 1
    while(r!=0)
    {
      a <- runif(1,0.5,0.6)
      u <- rbinom(P,1,a)
      u[which(u==0)] <- -1
      Vect_sign_dev <- u
      Dev <- Vect_sign_dev*Rvect
      MatDev <- matrix(rep(Beta+Dev,k-1),ncol=k-1)
      MatTest <- MatBeta-MatDev
      uTest <- apply(MatTest,2,ft)
      r <- length(which(uTest==0))
      if(r==0)
      {
        MatBeta <- cbind(MatBeta,MatDev[,1])
      }
      
    }
  }
  
  Vi <- runif(G,0.2,2)
  Ve <- runif(G,0.2,2)
  Vr <- runif(G,0.2,2)
  
  k <- 99
  while(k!=100)
  {
    ech <- sample(1:100,500,replace=T)  
    u_ech <- as.vector(table(ech))
    k <- length(u_ech)
  }
  n_ech <- rep(5,100)+u_ech
  
  Lw <- list(n_ech[1:33],n_ech[34:66],n_ech[67:100])
  To <- 10
  res <- SiMMMEI(MatBeta,Vi,Vr,Ve,Lw,To)
  
  Data <- res$Data
  MeData[[kk]] <- res
  formula <- YF~V2+V3+V4+V5+V6
  ResF <- EstimF(formula,Data,G)
  MeResF[[kk]] <- ResF
  
  uInit <- c(rep(1,33),rep(2,33),rep(3,34))
  
  
  Ind_R <- Ind_Regr(ResF$Rg,uInit)
  I1 <- c(I1,Ind_R[1])
  I2 <- c(I2,Ind_R[2])
  It <- c(It,ResF$Rs)
  D_Beta <- c(D_Beta,RMS_Beta(MatBeta[,unique(ResF$Rg)],ResF$theta$MatBeta))
  D_Vi <- c(D_Vi,RMS_Sigma(Vi[unique(ResF$Rg)],ResF$theta$Vi))
  D_Vr <- c(D_Vr,RMS_Sigma(Vr[unique(ResF$Rg)],ResF$theta$Vr))
  D_Ve <- c(D_Ve,RMS_Sigma(Ve[unique(ResF$Rg)],ResF$theta$Ve))
  print(kk)
}

mean(I1)
sd(I1)

mean(I2)
sd(I2)

mean(It)
sd(It)

mean(D_Beta)
sd(D_Beta)

mean(D_Vi)
sd(D_Vi)

mean(D_Ve)
sd(D_Ve)

mean(D_Vr)
sd(D_Vr)
