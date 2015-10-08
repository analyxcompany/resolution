#' cluster_resolution_RandomOrderFULL
#'
#' @description   cluster_resolution_RandomOrderFULL is enlargement of function cluster_resolution in the case where outcome based on random
#' order of nodes. Function repeats rep times algorithm and return all the obtained results (see the return section).
#' @param graph An igraph network or a data frame of three columns: source, target, and weights.
#' @param t The time-scale parameter of the process which uncovers community structures at different resolutions.
#' @param directed Logical. TRUE if the network is directed. Ignored if graph is an igraph object.
#' @param rep You can choose the number of repetitions. The default is 10.
#' @return  return list where we have (1) table with each result, (2) modularity for each outcome, (3) the best clustering,
#' (4) the the highest value of modularity which the best clustering has.
#' @examples
#' library(igraph)
#' g <- nexus.get("miserables")
#' cluster_resolution_RandomOrderFULL(g,directed=FALSE,t=1)
#' cluster_resolution_RandomOrderFULL(g,directed=FALSE,t=1,rep=20)

cluster_resolution_RandomOrderFULL <- function(graph, t = 1, directed=FALSE,rep=10)
{

  if(igraph::is.igraph(graph)){
    allVertex <- igraph::get.vertex.attribute(graph)$name
    g <- graph
  } else {
    graph <- as.data.frame(graph)
    allVertex <- unique(c(graph[,1],graph[,2]))
    if(length(which(graph[,3]==0))>0){
      graph <- graph[-which(graph[,3]==0),]}
    g <- igraph::graph.data.frame(graph, directed=directed)
  }

  A <- igraph::get.adjacency(g, type="both",
                             attr=names(edge.attributes(g)), edges=FALSE, names=TRUE,
                             sparse=FALSE)


  sampleorderL <- list()
  A_L <- list()
  for(i in 1:rep){
    sampleorderL[[i]] <- sample(rownames(A))
    A_L[[i]] <- A[sampleorderL[[i]],]
    A_L[[i]] <- A_L[[i]][,sampleorderL[[i]]]
  }
  sampleorderL <<- sampleorderL


  NodesGroupsT <- data.frame()
  for(i in 1:length(A_L)){
    CM <- as.matrix(A_L[[i]])

    #Initializing the table that informs about the community number to which a particular node is assigned
    NodesGroups <- data.frame(community=1:nrow(CM))
    rownames(NodesGroups) <- rownames(CM)

    # computing the matrix e^(t(B-I))k_j
    Adj <- Exp_matrix(CM,t)
    colnames(Adj) <- rownames(Adj)

    logic <- 1

    while (logic != 0 & nrow(CM) > 1 )
    {
      # names contains the names of communities which have not been "visited" yet in this iteration
      # a community is "visited" if it is merged with some other community or it was confirmed that no
      # increase in modularity  can be achieved by merging this community with any other
      names <- rownames(CM)

      # logic = 1 means that some increase of modularity has been achieved in the iteration
      logic <- 0

      while (length(names) > 0)
      {
        name1 <- names[1]
        max <- 0
        # analyze all the nodes in the neighborhood of name1
        neighborhood <- setdiff(rownames(CM)[which(CM[name1,]!=0)],name1)

        if(length(neighborhood) > 0)
        {
          delta <- Delta(CM,Adj,name1,neighborhood)
          if(length(neighborhood)==1){ names(delta) <- neighborhood}
          if(max(delta) > max){max <- max(delta); where<- names(delta)[which.max(delta)]; logic<-1}
        }


        # we have "visited" community name1 so we remove it from names
        names <- setdiff(names,name1)
        if(max >0)
        {
          # If there has been an increase in modularity
          # perform the merging that maximize this increase
          # modify the matrix CM and Adj
          # modify the NodesGroups table
          CM <- Erase(CM,name1,where)
          Adj <- Erase(Adj,name1,where)
          NodesGroups <- NewNodesGroup(name1,where,NodesGroups)

          #now we have also "visited" community "where" so we remove it from names
          names <- setdiff(names,where)
          if(length(CM) == 1) {CM <- as.matrix(CM)}
        }
      }
    }

    # modify the numbers denoting communities so that they are consecutive natural numbers starting from 1
    NodesGroups$community <- as.factor(NodesGroups$community)
    levels(NodesGroups$community) <- 1:length(levels(NodesGroups$community))
    NodesGroups$community <- as.numeric(NodesGroups$community)

    if(i==1){
      NodesGroupsT <- data.frame(name=rownames(NodesGroups),c1=NodesGroups)
      colnames(NodesGroupsT)[ncol(NodesGroupsT)] <- paste("c",i,sep="")
    } else {
      NodesGroups <- data.frame(name=rownames(NodesGroups),com=NodesGroups)
      NodesGroupsT <- merge(NodesGroupsT,NodesGroups,by="name",all=TRUE)
      colnames(NodesGroupsT)[ncol(NodesGroupsT)] <- paste("c",i,sep="")
    }
  }

  mod <- c()
  v <- V(g)$name
  for(i in 1:(ncol(NodesGroupsT)-1)){
    mod[i] <- igraph::modularity(g,NodesGroupsT[,i+1][match(v,NodesGroupsT[,1])],weights = E(g)$weight) }
  names(mod) <- colnames(NodesGroupsT)[-1]
  mod <<-mod


  if(length(setdiff(allVertex,NodesGroupsT$name))>0){
    temp <- data.frame(name=setdiff(allVertex,NodesGroupsT$name))
    NodesGroupsT <- merge(NodesGroupsT,temp,by="name",all=TRUE)
  }

  Results <- list()
  Results$communities <- NodesGroupsT
  Results$modularity <- mod
  Results$theBestPartition <- NodesGroupsT[,c(1,which.max(mod)+1)]
  Results$theHighestModularity <- mod[which.max(mod)]


  return (Results)
}
