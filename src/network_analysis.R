# Network analysis with random 200 nodes

# Install the required libraries
if(!"RCy3" %in% installed.packages()){
  install.packages("BiocManager")
  BiocManager::install("RCy3")
}

if(!"igraph" %in% installed.packages()){
  install.packages("BiocManager")
  BiocManager::install("igraph")
}

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
set.seed(135)
cmd.string = 'string protein query query="TP53" cutoff=0.9 species="Homo sapiens" limit=199'
commandsRun(cmd.string)

# Transfer network from Cytoscape to R after renaming it on Cytoscape
network <- createIgraphFromNetwork("TP53")

# Decide the layout
## A list of available layouts
getLayoutNames()

## Select the “force-directed” layout. To see properties for the given layout, use:
getLayoutPropertyNames("force-directed") 
layoutNetwork('force-directed defaultSpringLength=70 defaultSpringCoefficient=0.000003')

# The number of vertices and edges
igraph::vcount(network)
igraph::ecount(network)

# The names of vertices and edges
V(network)
E(network)

# Mean distance
mean_distance(network, directed = FALSE)
[1] 1.938408

# Diameter
diameter(network)

# Distances 
dist <- distances(network)

# Articulation point
## Removal of this point would destroy the whole network
articulation.points(network)
ENSP00000269305(TP53)

# Degree distribution
deg <- degree(network, mode="all")
hist(deg, breaks=1:vcount(network)-1, 
     ylim = range(pretty(c(0,table(deg)))),
     main="Histogram of Node Degree",
     xlab = "Degree")

# Betweennes centrality
bet <- betweenness(network)

bet[which.max(bet)]
ENSP00000269305(TP53)

# Closeness centrality
clo <- closeness(network, normalized = T)

clo[which.max(clo)]
ENSP00000269305(TP53)

# Save the session and export
full.path=paste(getwd(),'random_network_with_200_nodes',sep='/')
saveSession(full.path) #.cys

# Save image files with high resolution
full.path=paste(getwd(),'random_network_with_200_nodes',sep='/')
exportImage(full.path, 'PNG', zoom=500) #.png scaled by 500%
exportImage(full.path, 'PDF') #.pdf

# Track versions for records
cytoscapeVersionInfo()
sessionInfo()
