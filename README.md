# resolution
This is the R implementation of an algorithm to find communities in networks with resolution parameter based on the article ["Laplacian dynamics and Multiscale Modular Structure in Networks" R. Lambiotte et al.](http://arxiv.org/pdf/0812.1770.pdf)

The algorithm in the function **cluster_resolution** detects clusters using stability as an objective function to be optimised in order to find the best partition of network. The number of clusters typically decreases as the resolution parameter (t) grows, from a partition of one-node communities which are as many as nodes when t = 0 to a two-way partition as t grows.

Because of fact that the result of the algorithm depends on the order in which the nodes are considered, we have introduced the parameter RandomOrder. When RandomOrder is NULL the order of vertices come from the graph, if is FALSE vertices will be arranged in alphabetical order, and if is TRUE vertices will be arrange in random order. In the case of a choice random order you can set the number of repetitions (rep patameter) and then the best solution (which will have the highest value of modularity) among these repetitions will be returned. 

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

For this example we use the coappeareance network from Les Miserables, by Victor Hugo. Get more details about this data set with `igraph::nexus.info("miserables")`

```R
library(igraph) 
g <- nexus.get("miserables")
cluster_resolution(g)
cluster_resolution_RandomOrderFULL(g)
```

**data frame input**

For those not familiar with the 'igraph' package, is possible to calculate the algorithm directly from a data frame. This data frame should consist in three columns: 'from', 'to', and 'weights', indicating the corresponding nodes connections and the weights.

```R
data <- get.data.frame(g)
```

After you have your data in that format, the applications of the functions are equivalent to the previous ones, with one exception, the parameter directed indicating if the network is directed or not (directed = FALSE by default).
In this example, the miserables is an undirected graph so we don’t change it.

```R
cluster_resolution(data)
cluster_resolution_RandomOrderFULL(data)
```

**cluster_resolution**

In this function we can change the resolution parameter t, RandomOrder parameter (FALSE- alphabetical order, TRUE- random order, NULL- order from graph), number of repetition if RandomOrder is TRUE or directed parameter if the graph is directed.

```R
cluster_resolution(g)
cluster_resolution(g, t=0.5)
cluster_resolution(g,RandomOrder=NULL)
cluster_resolution(g,RandomOrder=FALSE)
cluster_resolution(g,RandomOrder=TRUE)
cluster_resolution(g,RandomOrder=TRUE,rep=10)
```

If we use igraph input function the output is a "communities" class object (as the cluster algorithms implemented in igraph package return). Please look the example below or see the http://igraph.org/r/doc/communities.html to find out what you can do.

If we use data.frame input we receive single column table with informations about community which has been found for each node. Rownames contains information about nodes.

```R
c <- cluster_resolution(g,directed=FALSE,t=1,RandomOrder=TRUE,rep=3)
c$membership  # A numeric vector, one value for each vertex, the id of its community.
c$memberships # It returns all the obtained results in matrix where columns corespond to the vertices and rows to the repetitions.
c$modularity  # Vector of modularity for each reperitions.
c$names       # Names od nodes.
c$vcount      # How many communities have been founded.
c$algorithm   # The name of the algorithm that was used to calculate the community structure
print(c)      # Prints a short summary.
membership(c) # The (best) membership vector, which had the highest value of modularity.
modularity(c) # The highest modularity value.
length(c)     # The number of communities.
sizes(c)      # Returns the community sizes, in the order of their ids.
algorithm(c)  # The name of the algorithm that was used to calculate the community structure.
crossing(c,g) # Returns a logical vector, with one value for each edge, ordered according to the edge ids. The value is TRUE iff the edge connects two different communities, according to the (best) membership vector, as returned by membership().

```

**cluster_resolution_RandomOrderFULL**

In this function we can change the resolution parameter t, the number of repetition (the default is 10) or directed parameter if the graph is directed. You can use it in both cases (igraph or data.frame input) and the output is four-element list containing: table with each outcome, modularity for each outcome, the best clustering (partition which has the highest value of modularity), the value of modularity fot the best clustering.

```R
cluster_resolution_RandomOrderFULL(g)
cluster_resolution_RandomOrderFULL(g,t=0.5)
cluster_resolution_RandomOrderFULL(g,rep=20)
```

