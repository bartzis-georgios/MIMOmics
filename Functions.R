#############################
######### Functions #########
#############################

### ipak function: installs and loads multiple R packages.
### check to see if packages are installed. Install them if they are not, then load them into the R session.

ipak <- function(pkg){
    new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
    if (length(new.pkg)) 
        install.packages(new.pkg, dependencies = TRUE)
    sapply(pkg, require, character.only = TRUE)
}

### "soft.thres" function for determining soft thresholding parameter.
### Takes as input a data matrix and computes a power gamma that sufficiently 
### supresses the magnitude of the random noise resulted from all p(p-1)/2
### correlations between Gaussian random variables bellow 1
soft.thres<-function(x){
    px<-dim(x)[2]
    Nx<-dim(x)[1]
    gamma<-1
    result<-(px*(px-1))/(2*((sqrt(Nx))^gamma))
    while (result>1){
           gamma<-gamma+1
           result<-(px*(px-1))/(2*((sqrt(Nx))^gamma))
           }
    return(gamma)
}


### "DensCentrHet" function for printing the Density, Centralization and Heterogeneity
### of each cluster in the network.Density denotes what is the average strength between the metabolites.
### High values of centralization (near 1) denote that there is at least one metabolite in 
### the module with more than average connections. High heterogeneity denotes that the weights are 
### distributed unevenly across the edges. The input on the function is the similarity matrix and 
### the vector denoting the color codes of each metabolite
DensCentrHet1<-function(x,y){
  quantities<-matrix(0,dim(table(y))+1,4)
  rownames(quantities)<-c(names(table(y)),"overall")
  colnames(quantities)<-c("nodes","Density","Centralization","Heterogeneity")
  if(sum(rownames(quantities)=="grey")==1){
    quantities<-quantities[-which(rownames(quantities)=="grey"),]
  }
  
  p1<-dim(x)[1]
  dens<-sum(x)/(p1*(p1-1))
  sbar<-mean(rowSums(x))
  smax<-max(rowSums(x))
  centr<-(smax-sbar)/p1
  heter<-sqrt(var(rowSums(x)))/sbar
  
  quantities[dim(quantities)[1],1]<-p1
  quantities[dim(quantities)[1],2]<-dens
  quantities[dim(quantities)[1],3]<-centr
  quantities[dim(quantities)[1],4]<-heter
  for (color in rownames(quantities)[-dim(quantities)[1]]){
    CBS<-x[which(y==color),which(y==color)]
    
    p1<-dim(CBS)[1]
    dens<-sum(CBS)/(p1*(p1-1))
    sbar<-mean(rowSums(CBS))
    smax<-max(rowSums(CBS))
    centr<-(smax-sbar)/p1
    heter<-sqrt(var(rowSums(CBS)))/sbar
    
    quantities[color,1]<-p1
    quantities[color,2]<-dens
    quantities[color,3]<-centr
    quantities[color,4]<-heter
  }
  quantities<-rbind(quantities[order(-quantities[1:(dim(quantities)[1]-1),2]),],quantities[dim(quantities)[1],])
  rownames(quantities)[dim(quantities)[1]]<-"overall"
  quantities<<-quantities
  return(signif(quantities,digits=3))
}


### "topCpercent" function takes as input the similarity matrix and c (=the percentage of edges to be kept)
### and returns the adjacency matrix with the top c% of the edges. It also stores the result in the object
### topCmatrix
topCpercent<-function(x,c){
    q <- quantile(upperTriangle(x, diag=FALSE),1-c);
    topCmatrix<<-(x>q)*1
    return(topCmatrix)
}

### "mostConnected" function gets the matrix x with the top 10% of the edges and gives a list 
### with the top #nr connected metabolites in decreasing order
mostConnected<-function(x,nr){
    list<-cbind(colnames(x),rowSums(x))
    return(print(cbind(list[,1][order(-as.numeric(list[,2]))],sort(as.numeric(list[,2]),decreasing=TRUE))[1:nr,],row.names=FALSE,quote=FALSE))
}


### "SoftConnectivity" function computes a vector of soft connectivities of nodes in weighted networks.
SoftConnectivity <- function(x){
   y<-(abs(cor(x,use="pairwise.complete.obs")))^soft.thres(x)
   diag(y)<-0
   y<-rowSums(y)
   return(y)
}


### "HardConnectivity" function computes a vector of soft connectivities of nodes in weighted networks.
HardConnectivity <- function(x){
   y<-(abs(cor(x,use="pairwise.complete.obs")))^soft.thres(x)
   diag(y)<-0
   y<-topCpercent(x=y,c=0.1)
   y<-rowSums(y)
   return(y)
}


### Function for simple network visualization
  visualize.network<-function(adjacency.matrix){
  ### Recover a vector with the pairwise edges
  edges=graph.adjacency(adjmatrix=adjacency.matrix, mode="undirected")
  ### Recover degree of each node
  degree=apply(adjacency.matrix, 2, sum)
  ### Find which are the nodes belonging in the top 50% with highest degree
  top50perc<-which(degree>=quantile(degree,prob=0.5))
  ### Find which are the nodes belonging in the top 25% with highest degree
  top25perc<-which(degree>=quantile(degree,prob=0.75))
  ### Find which are the nodes belonging in the top 5% with highest degree
  top5perc<-which(degree>=quantile(degree,prob=0.95))
  ### Color all nodes blue
  V(edges)$color<-"blue"
  ### Color top 50% connected with green
  V(edges)$color[top50perc]<-"green"
  ### Color top % connected with yellow
  V(edges)$color[top25perc]<-"yellow"
  ### Color top 5% connected with red
  V(edges)$color[top5perc]<-"red"
  net.plot<-plot(edges, vertex.size=4, vertex.frame.color="white",layout=layout.fruchterman.reingold, vertex.label=NA, edge.color=grey(0.5))
  return(net.plot)
}


### Function for checking degree distribution
degree.distribution<-function(adjacency.matrix){
  degree.plot<-plot(table(rowSums(adjacency.matrix)))
  return(degree.plot)
}

    
### Function for writing the edges in the Gephi format
write.edges<-function(matrix){
### Take the upper triangle of the  adjacency matrix
matrix<-upper.tri(matrix)*matrix
### Calculate how many edges you have
edges.num<-sum(upperTriangle(matrix,diag=FALSE)!=0)
### Make an empty matrix for the edges
edges<-matrix(0,edges.num,4)
edges<-as.data.frame(edges)
### In the 1st column, put the node from which the connection starts
edges[,1]<-as.numeric(which(matrix!=0,arr.ind = T)[,1])
edges[,1]<-colnames(matrix[,edges[,1]])
### In the 2nd column, put the node to which the connections end
edges[,2]<-as.numeric(which(matrix!=0,arr.ind = T)[,2])
edges[,2]<-colnames(matrix[,edges[,2]])
### Put the weights of the connections
edges[,3]<-upperTriangle(matrix)[upperTriangle(matrix) != 0 ]
### Indicate that the edges are indirect
edges[,4]<-"Undirected"
write.table(edges,"edges.csv",quote=FALSE,row.names=FALSE,col.names=c("Source","Target","Weight","Type"), sep=",")
}