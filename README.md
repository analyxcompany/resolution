# resolution
This is the R implementations of the algorithm finding communities in network with resolution parameter which has been created based on
paper ["Laplacian dynamics and Multiscale Modular Structure in Networks" R. Lambiotte et al.](http://arxiv.org/pdf/0812.1770.pdf)

Algorithm in function **cluster_resolution** detects clusters using stability as an objective function to be optimised in order to find the best partition of network. The number of clusters typically decreases as time grows, from a partition of one-node communities which are as many as nodes when t = 0 to a two-way partition as t → ∞.

Because of fact that the result of the algorithm depends on the order in which the nodes are considered, we have introduced parameter RandomOrder which if is NULL we receive the outcome based on order of vertices come from graph, if is FALSE vertices will be arrange in alphabetical order, if is TRUE vertices will be arrange in random order. In the case of a choice random order you can set the number of repetitions (rep patameter) and then will be returned the best solution (which will have the highest value of modularity) among these repetitions. 

In order to receive all outcomes from random orders it has been created **cluster_resolution_RandomOrderFULL** function which returns four-element list containing: table with each outcome, modularity for each outcome, the best clustering (partition which has the highest value of modularity), the value of modularity fot the best clustering.

## Installation
This package is not yet available in CRAN, so install it directly from Github with:

```R
# install.packages("devtools")
devtools::install_github("analyxcompany/resolution")
```

## Usage

After installation the package is loaded as usual with:

```R
library(resolution)
```

This implementation accepts as inputs an `igraph` object or a data frame.

**igraph input**

For this example I will use the coappeareance network from Les Miserables, by Victor Hugo. Get more details about this data set with `igraph::nexus.info("miserables")`

```R
library(igraph) 
g <- nexus.get("miserables")
cluster_resolution(g)
cluster_resolution_RandomOrderFULL(g)
```


**data frame input**

For those not familiar with the 'igraph' package, is possible to calculate the algorithm directly from a data frame. 
This data frame should consist in three columns: 'from', 'to', and 'weights', indicating the corresponding nodes connections and 
the weights.

```R
data <- get.data.frame(g)
```

After you have your data in that format, the applications of the functions are equivalent to the previous ones, 
with one exception, the parameter directed indicating if the network is directed or not (directed = FALSE by default).
In this example, the miserables is an undirected graph so we don’t change it.

```R
cluster_resolution(data)
cluster_resolution_RandomOrderFULL(data)
```

**cluster_resolution**

In this function we can change the resolution parameter t, RandomOrder patameter (FALSE- alphabetical order, TRUE- random order, NULL- order from graph), number of repetition if RandomOrder is TRUE or directed parameter if the graph is directed.

```R
cluster_resolution(g)
cluster_resolution(g, t=0.5)
cluster_resolution(g,RandomOrder=NULL)
cluster_resolution(g,RandomOrder=FALSE)
cluster_resolution(g,RandomOrder=TRUE)
cluster_resolution(g,RandomOrder=TRUE,rep=10)
```
The output is a data.frame where in column are communities which have been found for each node. Rownames contains information
about nodes.

**cluster_resolution_RandomOrderFULL**

In this function we can change the resolution parameter t, number of repetition (the default is 10) or directed parameter if the graph is directed.

```R
cluster_resolution_RandomOrderFULL(g)
cluster_resolution_RandomOrderFULL(g,t=0.5)
cluster_resolution_RandomOrderFULL(g,rep=20)
```
The output is four-element list containing: table with each outcome, modularity for each outcome, the best clustering (partition which has the highest value of modularity), the value of modularity fot the best clustering.

