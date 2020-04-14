
# generic methods for "DpeakData" class

setGeneric( "printEmpty",
    function( object, ... )
    standardGeneric("printEmpty")
)

setGeneric( "dpeakFit",
    function( object, ... )
    standardGeneric("dpeakFit")
)

setGeneric(
  "get_fragSet",
  function(x) standardGeneric("get_fragSet")
)

setGeneric(
  "get_PET",
  function(x) standardGeneric("get_PET")
)

setGeneric(
  "get_fragLenTable",
  function(x) standardGeneric("get_fragLenTable")
)

setGeneric(
  "get_Fratio",
  function(x) standardGeneric("get_Fratio")
)

setGeneric(
  "get_aveFragLen",
  function(x) standardGeneric("get_aveFragLen")
)

setGeneric(
  "get_stackedFragment",
  function(x) standardGeneric("get_stackedFragment")
)

setGeneric(
  "get_peakChr",
  function(x) standardGeneric("get_peakChr")
)

setGeneric(
  "get_peakStart",
  function(x) standardGeneric("get_peakStart")
)

setGeneric(
  "get_peakEnd",
  function(x) standardGeneric("get_peakEnd")
)

setGeneric(
  "get_emptyList",
  function(x) standardGeneric("get_emptyList")
)

# generic methods for "DpeakMotif" class

setGeneric(
  "get_motif",
  function(x) standardGeneric("get_motif")
)

setGeneric(
  "get_locMotif",
  function(x) standardGeneric("get_locMotif")
)

# generic methods for "DpeakFit" class

setGeneric( "export",
    function( object, ... )
    standardGeneric("export")
)

setGeneric(
  "get_fits",
  function(x) standardGeneric("get_fits")
)

setGeneric(
  "get_optFit",
  function(x) standardGeneric("get_optFit")
)

setGeneric(
  "get_optMu",
  function(x) standardGeneric("get_optMu")
)

setGeneric(
  "get_optPi",
  function(x) standardGeneric("get_optPi")
)

setGeneric(
  "get_optPi0",
  function(x) standardGeneric("get_optPi0")
)

setGeneric(
  "get_optGamma",
  function(x) standardGeneric("get_optGamma")
)

setGeneric(
  "get_optDelta",
  function(x) standardGeneric("get_optDelta")
)

setGeneric(
  "get_optSigma",
  function(x) standardGeneric("get_optSigma")
)

setGeneric(
  "get_bicVec",
  function(x) standardGeneric("get_bicVec")
)

setGeneric(
  "get_aicVec",
  function(x) standardGeneric("get_aicVec")
)

setGeneric(
  "get_estDeltaSigma",
  function(x) standardGeneric("get_estDeltaSigma")
)

setGeneric(
  "get_nTop",
  function(x) standardGeneric("get_nTop")
)

setGeneric(
  "get_lbDelta",
  function(x) standardGeneric("get_lbDelta")
)

setGeneric(
  "get_lbSigma",
  function(x) standardGeneric("get_lbSigma")
)

setGeneric(
  "get_psize",
  function(x) standardGeneric("get_psize")
)

setGeneric(
  "get_maxComp",
  function(x) standardGeneric("get_maxComp")
)

setGeneric(
  "get_pConst",
  function(x) standardGeneric("get_pConst")
)

setGeneric(
  "get_iterInit",
  function(x) standardGeneric("get_iterInit")
)

setGeneric(
  "get_iterMain",
  function(x) standardGeneric("get_iterMain")
)

setGeneric(
  "get_epsilon",
  function(x) standardGeneric("get_epsilon")
)
