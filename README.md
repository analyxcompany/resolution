# resolution
This is the R implementation of the algorithm finding communities in network with resolution parameter which has been created based on
paper ["Laplacian dynamics and Multiscale Modular Structure in Networks" R. Lambiotte et al.](http://arxiv.org/pdf/0812.1770.pdf)

Algorithm detects clusters using stability as an objective function to be optimised in order to find the best partition  of
network. The number of clusters typically decreases as time grows, from a partition of one-node communities which are
as many as nodes when t = 0 to a two-way partition as t → ∞.

## Installation
This package is not yet available in CRAN, so install it directly from Github with:

```R
# install.packages("devtools")
devtools::install_github("analyxcompany/resolution")
```

## Usage

After installation the package is loaded as usual with:

```R
library(ForceAtlas2)
```

This implementation accepts as inputs an `igraph` object or a data frame.

*igraph input**

For this example I will use the coappeareance network from Les Miserables, by Victor Hugo. Get more details about this data set with `igraph::nexus.info("miserables")`

```R
library(igraph) 
g <- nexus.get("miserables")
cluster_resolution(g)
cluster_resolution(g,t=0.5)
```

The default t parameter is 1 if you want you can change it (see second example above). 

**data frame input**
For those not familiar with the 'igraph' package, is possible to calculate the algorithm directly from a data frame. 
This data frame should consist in three columns: 'from', 'to', and 'weights', indicating the corresponding nodes connections and 
the weights.

```R
data <- get.data.frame(g)
```

After you have your data in that format, the application of the function is equivalent to the previous one, 
with one exception, the parameter directed indicating if the network is directed or not (directed = FALSE by default).
In this example, the miserables is an undirected graph so we don’t change it.

```R
cluster_resolution(data)
cluster_resolution(data,t=0.5)
```
**output**

The output is a data.frame where in column are communities which have been found for each node. Rownames contains information
about nodes.
