################################
# Constructor.
setClass('cellLinesObj',
         representation(synId = 'data.frame',  					# a data frame containing the synIds!
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
################################
