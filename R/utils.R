#' Erase
#' @description Function for the modification of the f.i. adjacency matrix resulting from merging communities OldCom and NewCom,
#' the row and column corresponding to OldCom are erased.
#' @param A Matrix which will be modified.
#' @param OldCom the community which the node have been.
#' @param NewCom the community which the node moved to.


Erase <- function(A,OldCom,NewCom)
{
  A[NewCom,] <- A[NewCom,]+ A[OldCom,]
  A[NewCom,NewCom] <- A[NewCom,NewCom] + A[OldCom,OldCom]
  A[,NewCom] <- A[NewCom,]
  return(A[-which(rownames(A)==OldCom),-which(colnames(A)==OldCom)])
}

#' NewNodesGroup
#' @description Updates table NodesGroups.
#' @param OldCom the community which the node have been.
#' @param NewCom the community which the node moved to.
#' @param NodesGroups table stores the number of the community to which a particular node is currently assigned.

NewNodesGroup <- function(OldCom,NewCom,NodesGroups)
{
  OldNumber <- NodesGroups[OldCom,"community"]
  NewNumber <- NodesGroups[NewCom,"community"]
  NodesGroups[which(NodesGroups$community==OldNumber),"community"] <- NewNumber

  return(NodesGroups)
}

#' Exp_matrix
#' @description Function for creating the matrix e^(t(B-I))k_j which is needed to computed stability.
#' @param A Adjacency matrix
#' @param p resolution parameter t

Exp_matrix <- function(A,p)
{

  v <- colSums(A)
  v[which(v==0)] <- min(v[v!=0])/1000
  B <- as.matrix(A) %*% diag(1/v)

  Adj <- as.matrix(expm::expm(p*(B-diag(1,nrow(B))),method = "Pade",order = 8))
  Adj <- Adj %*% diag(colSums(A))

}

#' Delta
#'  @description Computes value Delta(R_NL) for each neighbor at the same time and returned vector delta with these values.
#'
#' Delta(R_NL) evaluates the change of stability (R_NL) by removing Com1 from its community and
#' then by moving it into a neighbouring community. The node Com1 is then in cluster_resolution algotithm placed in the community for which
#' this gain is maximum, but only if this gain is positive.
#' @param A Adjacency matrix
#' @param Adj Matrix calxulated from pattern e^(t(B-I))k_j
#' @param Com1 node which is considered to change the community
#' @param Com2 neighbors of node Com2

Delta<- function(A,Adj,Com1,Com2)
{

  k_i_in <- Adj[Com1,Com2]
  m <- sum(A)/2
  if(length(Com2)==1){ S_tot <- sum(A[Com2,]) } else { S_tot <- rowSums(A[Com2,])}
  k_i <- sum(A[Com1,])

  return (k_i_in/(m) -(k_i *S_tot)/(2*(m^2)) )
}
