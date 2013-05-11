# Use cell lines

require(synapseClient)
require(affy)
synapseLogin()

source("/Users/fredcommo/Documents/MyProjects/virtualIC50/analysis/data_utils_2.R")
source("/Users/fredcommo/Documents/MyProjects/virtualIC50/analysis/celline_2_tcga_pipeline.R")

.getDataEntity <- function(synId){
  e <- loadEntity(synId)
  Data <- read.csv(file.path(e$cacheDir, e$files), header = TRUE, sep = ' ',
                   stringsAsFactors = FALSE)
  return(Data)
}

sanger <- loadEntity('syn1690438') # Charles' synId
names(sanger$objects$sanger_data)

#########################################

# Curate Sanger cells information
Info <- sanger$objects$sanger_data$sanger_info
Info <- cbind.data.frame(row.names = rownames(Info), cosmicId = Info[,2],
                         cancerType = Info[,3], origin = Info[,4])
Info <- Info[order(rownames(Info)),]

# Curate Sanger AUC and IC50
# ent <- loadEntity('syn1419911') # From the jira page at https://sagebionetworks.jira.com/wiki/display/METAGENOMICS/Cancer+Cell+Line+Data
AUC <- .getDataEntity('syn1836925')       # Last release, In sock 2013-05-09
AUC <- AUC[order(rownames(AUC)), order(colnames(AUC))]
IC50 <- .getDataEntity('syn1836934')      # Last release, In sock 2013-05-09
IC50 <- IC50[order(rownames(IC50)), order(colnames(IC50))]

cellNames <- toupper(gsub('\\W', '', rownames(Info)))
ICnames <- toupper(gsub('\\W', '', colnames(IC50)))
ICmatch <- lapply(1:nrow(Info), function(i){
  cellName <- cellNames[i]
  matchName <- NA
  cat(i, cellName, '\t')
  if(any(grep(paste0('^', cellName, '$'), ICnames)))
    matchName <- colnames(IC50)[grep(paste0('^', cellName, '$'), ICnames)]
  cat(matchName, '\n')
  return (data.frame(cellName = rownames(Info)[i], ICid = matchName))
})
ICmatch <- do.call(rbind, ICmatch)
dim(ICmatch); head(ICmatch, n = 20)

# curate Sanger drug information
drugs <- sanger$objects$sanger_data$sanger_drug_info
drugs <- drugs[order(drugs$Name),]
drugNames <- toupper(gsub('\\W', '', drugs$Name))
drugMatch <- lapply(1:length(drugNames), function(i){
  drugName <- drugs$Name[i]
  matchName <- NA
  cat(i, drugName, '\t')
  if(any(grepl(paste0('^',drugNames[i], '$'), rownames(IC50))))
    matchName <- rownames(IC50)[grepl(paste0('^',drugNames[i], '$'), rownames(IC50))]
  cat(matchName, '\n')
  return(data.frame(drugName = drugName, shortName = matchName))
})
drugMatch <- do.call(rbind, drugMatch)
dim(drugMatch); head(drugMatch)
drugs <- cbind.data.frame(drugs[,1:2], IC50id = drugMatch$shortName, drugs[,-c(1:2)])

# Curate GE and add GE ids to Info
sangerGE <- exprs(getSangerExpr_MetaGenomics()) # Curated by Brig, modified by Brian 2012-09-04, syn427896
sangerGE <- sangerGE[ ,order(colnames(sangerGE))]

GEnames <- toupper(gsub('_.*', '', colnames(sangerGE)))
cellMatchGE <- lapply(1:length(cellNames), function(i){
  cellName <- as.character(cellNames[i])
  matchName <- NA
  cat(i, cellName, '\t')
  if(any(grepl(paste0('^',cellName, '$'), GEnames)))
    matchName <- colnames(sangerGE)[grep(paste0('^',cellName, '$'), GEnames)]
  cat(matchName, '\n')
  return(data.frame(cellName = cellName, GEname = matchName))
})
cellMatchGE <- do.call(rbind, cellMatchGE)
dim(cellMatchGE); head(cellMatchGE)

# Curate CNV and add CNV ids to Info
sanger <- loadEntity('syn1690438')
sangerCNV <- sanger$objects$sanger_data$sanger_cnv  # syn1417761
colnames(sangerCNV)[12] <- '7860_KIDNEY' # wrong annotation: was annotated 786O_KIDNEY instead of 7860_KIDNEY
sangerCNV <- sangerCNV[,order(colnames(sangerCNV))]

CNVnames <- toupper(gsub('_.*', '', colnames(sangerCNV)))
cellMatchCNV <- lapply(1:length(cellNames), function(i){
  cellName <- as.character(cellNames[i])
  matchName <- NA
  cat(i, cellName, '\t')
  if(any(grepl(paste0('^',cellName, '$'), CNVnames)))
    matchName <- colnames(sangerCNV)[grep(paste0('^',cellName, '$'), CNVnames)]
  cat(matchName, '\n')
  return(data.frame(cellName = cellName, CNVname = matchName))
})
cellMatchCNV <- do.call(rbind, cellMatchCNV)
dim(cellMatchCNV); head(cellMatchCNV)

# Curate mutations
mut <- .getDataEntity('syn1836914')
mut <- mut[order(rownames(mut)), order(colnames(mut))]
mutNames <- toupper(gsub('_.*', '', colnames(mut)))
cellMatchMut <- lapply(1:length(cellNames), function(i){
  cellName <- as.character(cellNames[i])
  matchName <- NA
  cat(i, cellName, '\t')
  if(any(grepl(paste0('^',cellName, '$'), mutNames)))
    matchName <- colnames(mut)[grep(paste0('^',cellName, '$'), mutNames)]
  cat(matchName, '\n')
  return(data.frame(cellName = cellName, mutName = matchName))
})
cellMatchMut <- do.call(rbind, cellMatchMut)
dim(cellMatchMut); head(cellMatchMut)


# Add ICids, GE ids, CNV ids and mutation ids to Info
fullInfo <- cbind.data.frame(Info, drugTestId = ICmatch$ICid, GEid = cellMatchGE$GEname,
                             CNVid = cellMatchCNV$CNVname, MUTid = cellMatchMut$mutName)
synIds <- data.frame(Data = c('GeneExpr', 'CopyNumber', 'Mutations', 'IC50', 'AUC', 'DrugInfo', 'CellsInfo'),
                     synId = c('syn427896', 'syn1417761', 'syn1836914', 'syn1836934', 'syn1836925', NA, NA))

#########################################
# Build a Robject

# Constructor.
setClass('cellLinesObj',
         representation(synId = 'data.frame',						# a data frame containing the synIds!
                        exprSet = 'data.frame',					# a data.frame containing the expression set
                        cnvSet = 'data.frame', 					# a data.frame containing the CNV.
                        mutSet = 'data.frame',   				# a data.frame containing the mutations.
                        IC50 = 'data.frame',				    # a data.frame containing the IC50.
                        AUC = 'data.frame',					    # a data.frame containing the AUC50.
                        drugInfo = 'data.frame',  			# a data.frame containing information about drugs.
                        cellsInfo = 'data.frame')				# # a data.frame containing information about the cell lines.
      )

cellsObj <- function(synId, exprSet, cnvSet, mutSet, IC50, AUC, drugInfo, cellsInfo)
        {
        new('cellLinesObj', synId = synIds, exprSet = exprSet, cnvSet = cnvSet, mutSet = mutSet,
            IC50 = IC50, AUC = AUC, drugInfo = drugInfo, cellsInfo = cellsInfo)
        }

# Accessors
setGeneric("getInfo", function(object, arg = NULL,...) standardGeneric("getInfo"))
setGeneric("getExprs", function(object, arg = NULL,...) standardGeneric("getExprs"))
setGeneric("getCNV", function(object, arg = NULL,...) standardGeneric("getCNV"))
setGeneric("getMut", function(object, arg = NULL,...) standardGeneric("getMut"))
setGeneric("getIC50", function(object, arg = NULL,...) standardGeneric("getIC50"))
setGeneric("getAUC", function(object, arg = NULL,...) standardGeneric("getAUC"))
setGeneric("getDrugs", function(object, arg = NULL,...) standardGeneric("getDrugs"))
setGeneric("getCells", function(object, arg = NULL,...) standardGeneric("getCells"))

setMethod('getInfo', signature = 'cellLinesObj', function(object, arg = NULL,...){out = data.frame(object@synId); return(out)})
setMethod('getExprs', signature = 'cellLinesObj', function(object, arg = NULL,...){return(object@exprSet)})
setMethod('getCNV', signature = 'cellLinesObj', function(object, arg = NULL,...){return(object@cnvSet)})
setMethod('getMut', signature = 'cellLinesObj', function(object, arg = NULL,...){return(object@mutSet)})
setMethod('getIC50', signature = 'cellLinesObj', function(object, arg = NULL,...){return(object@IC50)})
setMethod('getAUC', signature = 'cellLinesObj', function(object, arg = NULL,...){return(object@AUC)})
setMethod('getDrugs', signature = 'cellLinesObj', function(object, arg = NULL,...){return(object@drugInfo)})
setMethod('getCells', signature = 'cellLinesObj', function(object, arg = NULL,...){return(object@cellsInfo)})

# ShowMethod
setGeneric("show", function(object, arg = NULL,...) standardGeneric("show"))
setMethod('show', signature = 'cellLinesObj',
          function(object){
            out <- rbind.data.frame(dim(getExprs(object)), dim(getCNV(object)),
                                    dim(getMut(object)), dim(getIC50(object)),
                                    dim(getAUC(object)), dim(getDrugs(object)),
                                    dim(getCells(object)))
            colnames(out) <- c('nRow', 'nCol')
            out <- cbind.data.frame(getInfo(object), out)
            cat('\nInstance of class', class(object), 'with', length(slotNames(object)), 'slot(s)', '\n\n')
            cat('Use getInfo(object) to get the synIds')
            cat('\n\nAccessor functions:
getExprs(object), getCNV(object), getMut(object), getIC50(object),
getAUC(object), getDrugs(object), getCells(object)\n\n')
            print(out)
            })

sangerCells <- cellsObj(synId = synIds, exprSet = as.data.frame(sangerGE), cnvSet = as.data.frame(sangerCNV),
                    mutSet = as.data.frame(mut), IC50 = as.data.frame(IC50), AUC = as.data.frame(AUC),
                    drugInfo = as.data.frame(drugs), cellsInfo = as.data.frame(Info))
# Store locally
setwd('/Users/fredcommo/Documents/MyProjects/CellLines/SangerCells/')
save(sangerCells, file = 'sangerCells.RData')


#############################################################
#
# Doesn't work yet!
#
#############################################################

# Store in synapse
# demo of how to create a file then use the uploaded file in a wiki
projectId <- "syn1855795"

# Store RData manually
myPath <- paste0(getwd(),'/')
myData <- Data(list(name = "sangerCells.RData", parentId = projectId))
myData <- synStore(myData) # then go to synapse and upload the .RData manually

# Store RObject constructor and accessors
myCode <- Code(list(name = "cellLinesObjClass.R", parentId = projectId))
myCode <- addFile(myCode, paste0(myPath, 'cellLinesObjClass.R'))
myCode <- synStore(myCode)

# Store Rcode
myCode <- Code(list(name = "buildSangerObj.R", parentId = projectId))
myCode <- addFile(myCode, paste0(myPath, 'buildSangerObj.R'))
myCode <- synStore(myCode)

# Add a provenance
myData <- getEntity('syn1855799')
used(myData) <- list(list(entity = myCode, wasExecuted=T),
              list(entity = getEntity('syn427896'), wasExecuted=F),
              list(entity = getEntity('syn1417761'), wasExecuted=F),
                   list(entity = getEntity('syn1836914'), wasExecuted=F),
                   list(entity = getEntity('syn1836934'), wasExecuted=F),
                   list(entity = getEntity('syn1836925'), wasExecuted=F))
myData <- synStore(myData)



