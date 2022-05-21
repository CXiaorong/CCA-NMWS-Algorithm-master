#CCA-NMWS
#1.First select the function(1)-function(14) function and run it.
#2.load data
#3.Initialization parameters.
#4.Iterative loop


#function(1): Calculate the actual number of edges in the network
cor_sum <- function(Like,x){  
  k <- length(x)
  result <- 0
  if(k!=1){
    for(i in 1:(k-1)){
      for(j in (i+1):k){
        temp <- as.numeric(Like[x[i],x[j]])
        result <- result+abs(temp)
      }
    }
  }
  else{
    result <- 0
  }
  return(result)
}


#function(2): fitness function--NMWS model
fitness <- function(A,Like,x){
  temp <- matrix(0,1,n)
  temp[x] <- 1
  y <- temp
  index <- which(y==1) 
  m1=nrow(A)
  a <- as.matrix(A[,index])
  a_indexsum <- rowSums(a)
  a_indexsum2=a_indexsum[which(a_indexsum>1)] 
  a_indexsum3=a_indexsum[which(a_indexsum>0)] 
  
  a_colsum2= colSums(a)
  a_colsum <- colSums(a)/m1 
  l1=length(a_indexsum3)
  l2=length(a_indexsum2)
  
  f_s=matrix(0,1,k+4)
  f_s[1,1:k]=index
  if(l1==0){ 
    f_s[1,k+2]=0
  }else{
    f_s[1,k+2]=-sum(a_indexsum2*(a_indexsum2-1))/(k*(k-1)*l1)
  }
  second_like=cor_sum(Like,index)/(k*k)
  f_s[1,k+1]=length(a_indexsum3)/m1
  f_s[1,k+3]=second_like
  f_s[1,k+4]=f_s[1,k+1]+f_s[1,k+2]+f_s[1,k+3]
  return(f_s)
}


#function(3): Compute the selection probability function.
select_order_fitness2 <- function(fit_vector){
  n <- length(fit_vector) 
  p <- matrix(0,1,n)
  s=sum(fit_vector)
  p=fit_vector/s
  p_cumsum <- cumsum(p)
  random_data <- runif(1)
  temp <- which(p_cumsum>=random_data)
  index1 <- temp[1]
  return(index1)
}

#function(4): crossover function 
crossover3 <- function(parent1,parent2,n){ 
  temp22=parent1
  temp <- matrix(0,1,n) 
  temp[parent1] <- 1
  parent1 <- temp
  temp <- matrix(0,1,n)
  temp[parent2] <- 1
  parent2 <- temp
  newpop <- matrix(0,1,n)
  index <- which((parent1+parent2)==2)  #Take the intersection
  newpop[index] <- 1
  parent1[index] <- 0
  parent2[index] <- 0
  temp <- which(parent1+parent2==1) #take union
  index <- sample(1:sum(parent1+parent2==1),sum(parent1+parent2)) #Shuffle the sequence numbers of the union genes of 2 individuals
  newpop[temp[index[1:(k-sum(newpop))]]]=1
  newpop <- which(newpop==1)#produce a new offspring
  
  beixuan<-c(1:n) 
  beixuan <- beixuan[!beixuan%in%newpop] 
  parent3=sample(beixuan,k)
  ptemp=matrix(0,k+1,k+5)
  ptemp[1,1:(k+4)]=fitness(SNVdata,second_like_data,newpop)
  index2=sample(1:k,1)
  for (ii in 1:k) {
    pa=newpop[1:k]
    pa[index2]=parent3[ii]
    ptemp[(ii+1),1:(k+4)]=fitness(SNVdata,second_like_data,pa)
  }
  I <- order(ptemp[,k+4],decreasing = 'T')
  ptemp2=ptemp[I,]
  return(ptemp2)
}

#function(5): Mutation function
mutation_SA <- function(SNVdata,second_like_data,popTemp,n,N){
  popTemp2=matrix(0,2,k+5)
  popTemp2[1,]=popTemp
  pop_i=popTemp[1:k]
  pop_j=mutation(pop_i,n) #eg:基因n=1,2,3 ->n=1,2,6
  popTemp2[2,1:(k+4)]=fitness(SNVdata,second_like_data,pop_j)
  if (popTemp2[2,(k+4)]>=popTemp2[1,(k+4)])
    pop_i <- pop_j
  return(pop_i)
}

#function(6): eg: gene n=1,2,3 ->n=1,2,6
mutation <- function(x,n){
  temp <- matrix(0,1,n)
  temp[x] <- 1
  x <- temp
  k <- sum(x)
  index1 <- round(runif(1,min = 1,max = n))#
  while(x[index1]==1){
    index1 <- round(runif(1,min = 1,max = n))
  }
  x_nonzero <- which(x==1)
  index2 <- x_nonzero[round(runif(1,min = 1,max = k))]
  m <- x  
  m[index1]=1
  m[index2]=0
  m <- which(m==1)
  return(m)
}

#function(7): select function
select3 <- function(pop3){
  I <- order(pop3[,k+4],decreasing = T) 
  pop3 <- pop3[I,]
  p <- vector(length = popsize)
  pop_next <- matrix(0, popsize,k+5)
  pop_next[1,] <- pop3[1,] 
  ss=sum(as.numeric(pop3[,k+4]))
  p=as.numeric(pop3[,k+4])/ss
  for(i in 2:popsize){  
    random_data=runif(1)
    p_cumsum <- cumsum(p)
    temp <- which(p_cumsum>=random_data)
    index1 <- temp[1]
    pop_next[i,] <- pop3[index1,]
  }
  return(pop_next)
}



#function(8): competition function
compareCGA3 <- function(p1,p2,k) { 
  n=nrow(p1)
  m=ncol(p1)
  n2=nrow(p2)
  childPop=matrix(0,n+1,m+n2+1)
  childPop[1:n,1:m]=p1
  ScoreE=matrix(0,1,n2) 
  for (j in 1:n2) { 
    childPop[n+1,1:m]=p2[j,]
    popTmp=childPop
    I <- order(popTmp[,k+5]) 
    ScoreE[j]=which(I ==(n+1)) #Ranking of opponent sets
  }
  return(ScoreE)
}

#function(9): Compare the magnitude of two numbers
compare <- function(a,b) {
  if(a<=b)
    return(1) 
  else
    return(2)
}


#function(10): GA2 function
GA2 <- function(pop){  
  I <- order(pop[,k+4],decreasing = T) 
  pop <- pop[I,]
  p <- vector(length = popsize)
  pop_next <- matrix(0, popsize,k+5)
  pop_next[1,] <- pop[1,] 
  fit_vector <-as.numeric(pop[,k+4]) 
  n=ncol(SNVdata)
  #The top 50% of the population is obtained through crossover
  f1 <- floor(popsize*(1-pm1))
  for(i in 2:popsize){ 
    index <- select_order_fitness2(fit_vector) 
    index1 <- index[1]
    index2 <- index[2]
    cro_pop <- crossover3(pop[index1,(1:k)],pop[index2,1:k],n)
    pop_next[i,1:(k+4)]=cro_pop[1,1:(k+4)]
  }
  #50% of individuals are acquired through mutation
  I <- order(pop_next[,k+4],decreasing = T) 
  pop_next <- pop_next[I,]
  for(i in f1:popsize){ 
      fit_vector <-as.numeric(pop_next[,k+4])
      index <- select_order_fitness2(fit_vector) 
      index1 <- index[1] 
      if(index1==1){
        index1=index1+1
      }
      n4 <- ncol(SNVdata)
      new_ind=mutation_SA(SNVdata,second_like_data,pop[index1,],n4,1)
      pop_next[i,1:(k+4)]=fitness(SNVdata,second_like_data,new_ind)
    # }
  }
  return(pop_next) 
}


#function(11): GA4 function
GA4 <- function(pop){  
  I <- order(pop[,k+4],decreasing = T) 
  pop <- pop[I,]
  p <- vector(length = popsize)
  pop_next <- matrix(0, popsize,k+5)
  
  #The top 20% of the population is obtained through crossover
  pop_next[1,] <- pop[1,] 
  fit_vector <-as.numeric(pop[,k+4]) 
  n=ncol(SNVdata)
  f1 <- floor(popsize*(1-pm2))
  for(i in 2:popsize){ 
    index <- select_order_fitness2(fit_vector) #选择适应度
    index1 <- index[1]
    index2 <- index[2]
    cro_pop <- crossover3(pop[index1,(1:k)],pop[index2,1:k],n)

    pop_next[i,1:(k+4)]=cro_pop[1,1:(k+4)]
  }

  #80% of individuals are acquired through mutation
  I <- order(pop_next[,k+4],decreasing = T) 
  pop_next <- pop_next[I,]
  for(i in f1:popsize){ 
      fit_vector <-as.numeric(pop_next[,k+4])
      index <- select_order_fitness2(fit_vector) 
      index1 <- index[1] 
      if(index1==1){
        index1=index1+1
      }
      n4 <- ncol(SNVdata)
      new_ind=mutation_SA(SNVdata,second_like_data,pop[index1,],n4,1)
      pop_next[i,1:(k+4)]=fitness(SNVdata,second_like_data,new_ind)
  }
  return(pop_next) 
}

#function(12): Cooperative pool evolution process.
GA5 <- function(pop){  

  I <- order(pop[,k+4],decreasing = T) 
  pop <- pop[I,]
  p <- vector(length = popsize)
  pop_next <- matrix(-1000, popsize,k+5)
  #The top 20% of the population is obtained through crossover
  fit_vector <-as.numeric(pop[,k+4]) 
  n=ncol(SNVdata)
  f1 <- nrow(pop)
  pop_next[1,] <- pop[1,] #keep the best individual
  for(i in 2:f1){ 
    index <- select_order_fitness2(fit_vector) 
    index1 <- index[1]
    index2 <- index[2]
    cro_pop <- crossover3(pop[index1,(1:k)],pop[index2,1:k],n)
    pop_next[i,1:(k+4)]=cro_pop[1,1:(k+4)]
  }
  #The rest is obtained by mutation
  fit_vector <-as.numeric(pop_next[,k+4])
  index <- select_order_fitness2(fit_vector) 
  index1 <- index[1] 
  Kgene=pop_next[index1,1:k]
  beixuan<-c(1:n) 
  beixuan <- beixuan[!beixuan%in%Kgene] 
  popTemp2=matrix(-1000,2,k+5)
  popTemp2[1,]=pop_next[index1,]
  temp=sample(1:k,1)
  for(i in 1:(n-k)){
    Kgene[temp] <- beixuan[i]
    f2 <- fitness(SNVdata,second_like_data,Kgene)
    popTemp2[2,1:(k+4)]=f2
    if(popTemp2[2,k+4]>popTemp2[1,k+4]){
      popTemp2[1,]=popTemp2[2,]
      break
    }
  }
  pop_next[index1,]=popTemp2[1,]
  return(pop_next)
}

#function(13): Initialize the population function
initPop<-function(SNVdata,second_like_data,k,n,popsize){
  pop <- matrix(-1000, popsize,k+5)
  for(i in 1:(popsize)){
    x <- 1:n
    temp <- sample(x,n) #Shuffle the order of n genes, and only take the first 1:k each time
    pop[i,1:(k+4)]=fitness(SNVdata,second_like_data,temp[1:k])
  }
  return(pop) 
}

#function(14): significance test function
significance <- function(A,Like,subset_M){
  m <- nrow(A)
  W <- matrix(0,1,1000)
  n <- length(subset_M)
  for(j in 1:1000){
    A_temp <- A
    A_temp[,subset_M] <- 0
    for(i in 1:n){
      temp <- sum(A[,subset_M[i]])
      index <- round(runif(temp,min = 1,max = m))
      A_temp[index,subset_M[i]] <- 1
    }
    W1=fitness(A_temp,Like,subset_M)
    W[j]=W1[1,k+4]
  }
  W2=fitness(A,Like,subset_M)
  p <- sum(W>=W2[1,k+4])/1000
  return(p)
}


############################################################################################
#1.load data
#GBM
SNV_data<-read.csv('E:/sourceFiles/data/GBM/SNVdata_440.csv')
second_like_data<-read.csv('E:/sourceFiles/data/GBM/network_440.csv')
rownames(second_like_data)=second_like_data[,1]
second_like_data=second_like_data[,-1]
second_like_data=as.matrix(second_like_data)

#OVCA
# SNV_data<-read.csv('E:/sourceFiles/data/OVCA/SNVdata_2547.csv')
# second_like_data<-read.csv('E:/sourceFiles/data/OVCA/network_2547.csv')
# rownames(second_like_data)=second_like_data[,1]
# second_like_data=second_like_data[,-1]
# second_like_data=as.matrix(second_like_data)


#THCA
# SNV_data<-read.csv('E:/sourceFiles/data/THCA/SNVdata_3420.csv')
# second_like_data<-read.csv('E:/sourceFiles/data/THCA/network_3420.csv')
# rownames(second_like_data)=second_like_data[,1]
# second_like_data=second_like_data[,-1]
# second_like_data=as.matrix(second_like_data)

#2.Data preprocessing.
rownames(SNV_data)<-SNV_data[,1]
SNV_data<-SNV_data[,-1]
snv_colsum<-colSums(SNV_data) 
SNVdata<-SNV_data[,snv_colsum>1]


#3.Initialization parameters
n <- ncol(SNVdata)
m <- nrow(SNVdata)
geneName <- colnames(SNVdata)
iteration <- 1000    #number of iterations
pm1=0.5
pm2=0.8
cta=0.2
k <-5   #k
popsize <- floor(log2(n^k))*2 #population size


#4. Iterative loop, v is the number of executions of the algorithm
for(v in 1:10){
  str=paste("Execute the ", v,"th time", sep = "")
  print(str)
  # Initialize the population
  result_pop=initPop(SNVdata,second_like_data,k,n,2*popsize)
  P1=result_pop[1:popsize,]
  P2=result_pop[(popsize+1):(2*popsize),]
  x <- 1:popsize
  j=1
  R=1
  t=0
  I=order(P1[,k+4],decreasing = 'T')
  Pbest=P1[I,]
  source_pop=matrix(-1000,popsize,k+5)
  betR=floor((1/k)*100)
  if(betR<20){
    betR=20
  }
  #Determine if competition is required.
  while(j<=iteration & R<=10) { 
    
    best_ind=Pbest[1,1:k]
    temp_pop <- as.matrix(P1[!duplicated(P1[,1:k]),,drop=T])
    temp_pop2 <- as.matrix(P2[!duplicated(P2[,1:k]),,drop=T])
    R1=nrow(temp_pop)
    R2=nrow(temp_pop2)
    #Determine whether to converge.
    if(nrow(temp_pop)*ncol(temp_pop)==(k+5)||nrow(temp_pop2)*ncol(temp_pop2)==(k+5)){
      break
    }
    
    if((j%%10==0)&(R1<=popsize/2 || R2<=popsize/2)){

      #competition between populations
      I=order(P1[,k+4],decreasing = 'T')
      P1=P1[I,]
      I=order(P2[,k+4],decreasing = 'T')
      P2=P2[I,]
      source_pop[1:(popsize/2),]=P1[1:(popsize/2),]
      source_pop[(popsize/2+1):popsize,]=P2[1:(popsize/2),]
      source_pop=GA5(source_pop)
      #Select opponent set individuals
      E1=base::unique(P1) 
      E2=base::unique(P2) 
      if(nrow(E2)>=popsize/5){
        l1=popsize*cta
      }else{
        l1=nrow(E2)
      }
      if(nrow(E1)>=popsize/5){
        l2=popsize*cta
      }else{
        l2=nrow(E1)
      }
      P1C=compareCGA3(P1,E2[1:l1,],k) 
      P2C=compareCGA3(P2,E1[1:l2,],k) 
      
      #Calculate the population after source_pop, mix the winning population with the resource population, and take the first popsize.
      temp=matrix(0,2*popsize,k+5)
      temp[1:popsize,]=source_pop 
      if(sum(P1C)<sum(P2C)){ 
        temp[(popsize+1):(2*popsize),]=P1
        P1=temp[1:popsize,]
      }else{
        temp[(popsize+1):(2*popsize),]=P2
        P2=temp[1:popsize,]
      }
      
      I=order(P1[,k+4],decreasing = 'T')
      Pbest=P1[I,]
      ll=length(which(best_ind==Pbest[1,1:k])==TRUE) 
      if(ll==k){
        R=R+1
      }else{
        R=0
      }
    }

    P1=GA2(P1)
    P2=GA4(P2)
    P1 <- select3(P1)
    P2 <- select3(P2)
    j=j+1
  }
  
  pop=matrix(0,2*popsize,k+5)
  pop[1:popsize,]=P1
  pop[(popsize+1):(2*popsize),]=P2
  I <- order(pop[,(k+4)],decreasing = T) 
  pop <- pop[I,]
  maxpop <- pop
  p=significance(SNVdata,second_like_data,maxpop[1,1:k])
  
  maxpop[,1:k] <- apply(pop[,1:k,drop=T],2,function(x) {geneName[x]})

  #output the optimal solution.
  print(maxpop[1,1:k])
  print(maxpop[1,(k+4)])
  v=v+1
}
