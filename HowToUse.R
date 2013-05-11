
################
# How to use it:
################

# Locally
setwd('/Users/fredcommo/Documents/MyProjects/CellLines/SangerCells/')
source('./cellLinesObjClass.R')
load('sangerCells.RData')
sangerCells

# Using synapse
require(synapseClient)
synapseLogin()

e <- loadEntity('syn1855802')
source(file.path(e$cacheDir, e$files))

e <- loadEntity('syn1855799')
load(file.path(e$cacheDir, e$files))
sangerCells

ge <- getExprs(sangerCells)
ge[1:10, 1:5]

cellList <- getCells(sangerCells)
head(cellList)