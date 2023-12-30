# Network analysis with random 200 nodes

# Install the required libraries
if(!"RCy3" %in% installed.packages()){
  install.packages("BiocManager")
  BiocManager::install("RCy3")
}

install.packages("igraph")


# Prerequisites
## In addition to these packages (RCy3, igraph), you will need:
## Cytoscape, which can be downloaded from http://www.cytoscape.org/download.php.

## There is also need for the STRING app to access the STRING database from within Cytoscape: * Install the STRING app from https://apps.cytoscape.org/apps/stringapp
## available in Cytoscape 3.7.0 and above; installApp('STRINGapp')  

# Load the libraries
library(igraph)
library(RCy3)

# Make sure to launch Cytoscape, and now connect R to Cytoscape
cytoscapePing()

# Create data 
cmd.string = 'string protein query query="TP53" cutoff=0.9 species="Homo sapiens" limit=199'
commandsRun(cmd.string)

# Transfer network from Cytoscape to R after renaming it on Cytoscape
network <- createIgraphFromNetwork()

# Decide the layout
## A list of available layouts
getLayoutNames()

## Select the “force-directed” layout
layoutNetwork('force-directed')

# The number of vertices and edges
igraph::vcount(network)
igraph::ecount(network)

# The names of vertices and edges
V(network)
E(network)

# Diameter
diameter(network)

# Distances 
dist <- distances(network)

# Articulation point
articulation.points(network)

# Degree distribution
deg <- degree(network, mode="all")
hist(deg, breaks=1:vcount(network)-1, 
     ylim = range(pretty(c(0,table(deg)))),
     main="Histogram of Node Degree",
     xlab = "Degree")

# Betweennes centrality
bet <- betweenness(network)
bet[which.max(bet)]

# Closeness centrality
clo <- closeness(network, normalized = T)
clo[which.max(clo)]

# Save the session and export
full.path=paste(getwd(),'out/random_network_with_200_nodes',sep='/')
saveSession(full.path) #.cys

# Save image files with high resolution
full.path=paste(getwd(),'out/random_network_with_200_nodes',sep='/')
exportImage(full.path, 'PNG', zoom=500) #.png scaled by 500%
exportImage(full.path, 'PDF') #.pdf

# Track versions for records
cytoscapeVersionInfo()
sessionInfo()
