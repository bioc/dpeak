
# generic methods for "DpeakFit" class

setMethod(
  f = "get_fits",
  signature = "DpeakFit",
  definition = function(x) x@fits
)

setMethod(
  f = "get_optFit",
  signature = "DpeakFit",
  definition = function(x) x@optFit
)

setMethod(
  f = "get_optMu",
  signature = "DpeakFit",
  definition = function(x) x@optMu
)

setMethod(
  f = "get_optPi",
  signature = "DpeakFit",
  definition = function(x) x@optPi
)

setMethod(
  f = "get_optPi0",
  signature = "DpeakFit",
  definition = function(x) x@optPi0
)

setMethod(
  f = "get_optGamma",
  signature = "DpeakFit",
  definition = function(x) x@optGamma
)

setMethod(
  f = "get_optDelta",
  signature = "DpeakFit",
  definition = function(x) x@optDelta
)

setMethod(
  f = "get_optSigma",
  signature = "DpeakFit",
  definition = function(x) x@optSigma
)

setMethod(
  f = "get_bicVec",
  signature = "DpeakFit",
  definition = function(x) x@bicVec
)

setMethod(
  f = "get_aicVec",
  signature = "DpeakFit",
  definition = function(x) x@aicVec
)

setMethod(
  f = "get_fragSet",
  signature = "DpeakFit",
  definition = function(x) x@fragSet
)

setMethod(
  f = "get_PET",
  signature = "DpeakFit",
  definition = function(x) x@PET
)

setMethod(
  f = "get_fragLenTable",
  signature = "DpeakFit",
  definition = function(x) x@fragLenTable
)

setMethod(
  f = "get_Fratio",
  signature = "DpeakFit",
  definition = function(x) x@Fratio
)

setMethod(
  f = "get_aveFragLen",
  signature = "DpeakFit",
  definition = function(x) x@aveFragLen
)

setMethod(
  f = "get_stackedFragment",
  signature = "DpeakFit",
  definition = function(x) x@stackedFragment
)

setMethod(
  f = "get_peakChr",
  signature = "DpeakFit",
  definition = function(x) x@peakChr
)

setMethod(
  f = "get_peakStart",
  signature = "DpeakFit",
  definition = function(x) x@peakStart
)

setMethod(
  f = "get_peakEnd",
  signature = "DpeakFit",
  definition = function(x) x@peakEnd
)

setMethod(
  f = "get_estDeltaSigma",
  signature = "DpeakFit",
  definition = function(x) x@estDeltaSigma
)

setMethod(
  f = "get_nTop",
  signature = "DpeakFit",
  definition = function(x) x@nTop
)

setMethod(
  f = "get_lbDelta",
  signature = "DpeakFit",
  definition = function(x) x@lbDelta
)

setMethod(
  f = "get_lbSigma",
  signature = "DpeakFit",
  definition = function(x) x@lbSigma
)

setMethod(
  f = "get_psize",
  signature = "DpeakFit",
  definition = function(x) x@psize
)

setMethod(
  f = "get_maxComp",
  signature = "DpeakFit",
  definition = function(x) x@maxComp
)

setMethod(
  f = "get_pConst",
  signature = "DpeakFit",
  definition = function(x) x@pConst
)

setMethod(
  f = "get_iterInit",
  signature = "DpeakFit",
  definition = function(x) x@iterInit
)

setMethod(
  f = "get_iterMain",
  signature = "DpeakFit",
  definition = function(x) x@iterMain
)

setMethod(
  f = "get_epsilon",
  signature = "DpeakFit",
  definition = function(x) x@epsilon
)

setMethod(
    f="show",
    signature="DpeakFit",
    definition=function( object ) {
        # extract objects

        optMu <- get_optMu(object)
        optPi0 <- get_optPi0(object)
        maxComp <- get_maxComp(object)

        # summary

        med_opt_mu <- median( unlist( lapply( optMu,
            function(x) ifelse( any(!is.na(x)), length(x), 0 )
        ) ) )
        pi0_vec <- unlist( optPi0 )
        med_explain <- round( 100 * median( 1 - pi0_vec, na.rm=TRUE ) ) / 100

        cat( "------------------------------------------------------------\n" )
        cat( "Summary: Dpeak model fitting (class: DpeakFit)\n" )
        cat( "------------------------------------------------------------\n" )
        cat( "- Maximum possible number of binding events in each peak: ",maxComp,"\n", sep="" )
        cat( "- Median number of binding events in each peak: ",med_opt_mu,"\n", sep="" )
        cat( "- Median explanation ratio: ",med_explain," %\n", sep="" )
        cat( "------------------------------------------------------------\n" )
    }
)

setMethod(
    f="print",
    signature="DpeakFit",
    definition=function( x ) {
        warning( "'print' method for 'DpeakFit' class is not supported yet." )
    }
)

setMethod(
    f="plot",
    signature=c("DpeakFit","missing"),
    definition=function( x, y, filename=NULL, plotType="fit",
        strand=FALSE, extension=1, smoothing=FALSE,
        threshold=1000, nsimul=10000, seed=12345, nCore=1, ... ) {

        pdf( filename )

        if ( identical(plotType,"fit") ) {

            # fitting results

            .plotFit( object=x, threshold=threshold,
                strand=strand, extension=extension, smoothing=smoothing )
        } else if ( identical(plotType,"GOF") ) {

            .plotGOF( object=x, nsimul=nsimul, seed=seed,
                extension=extension, smoothing=smoothing, nCore=nCore )
        } else if ( identical(plotType,"BIC") ) {
            # plot of BIC curve

            .plotBIC( object=x )
        } else {
            stop( "Error: Inappropriate 'plotType'. Please choose either 'fit', 'GOF', or 'BIC'." )
        }

        dev.off()
    }
)

setMethod(
    f="exportPeakList",
    signature="DpeakFit",
    definition=function( object, type=NA, filename=NA, ... ) {
        # error treatment: check invalid type

        if ( is.na(type) )
        {
            message( "Info: 'type' is not specified by the user." )
            message( "Info: 'type' is specified as 'bed' instead." )
            type <- "bed"
        }

        allType <- c("txt","gff","bed")
        invalidType <- TRUE
        for ( i in seq_len(length(type)) )
        {
            if ( length(which(!is.na(match(type[i],allType)))) == 0 )
            {
                invalidType <- FALSE
            }
        }
        if ( !invalidType )
        {
            message( "Info: 'type' incorrect." )
            message( "Info: 'type' is specified as 'bed' instead." )
            type <- "bed"
        }

        # error treatment: 'filename' not specified

        if ( is.na(filename) )
        {
            message( "Info: 'filename' is not specified by the user." )
            message( "Info: 'filename' is specified as 'bindingSiteList' instead." )
            filename <- paste("bindingSiteList.",type,sep="")
        }

        # extract objects

        optMu <- get_optMu(object)
        optPi <- get_optPi(object)
        peakChr <- get_peakChr(object)
        peakStart <- get_peakStart(object)
        peakEnd <- get_peakEnd(object)
        psize <- get_psize(object)

        # process fitting results for exporting

        peakName <- paste( peakChr, ":", peakStart, "-", peakEnd, sep="" )
        chrVec <- muVec <- strengthVec <- nameVec <- c()
        for ( i in seq_len(length(optMu)) ) {
            # error treatment: skip peaks with no fragments
            if ( is.na(get_fragSet(object)[[i]][1,1]) ) {
                next;
            }

            # stack data

            chrVec <- c( chrVec, rep( peakChr[i], length(optMu[[i]]) ) )
            muVec <- c( muVec, optMu[[i]] )
            strengthVec <- c( strengthVec, nrow(get_fragSet(object)[[i]]) * optPi[[i]] )
            nameVec <- c( nameVec, rep( peakName[i], length(optMu[[i]]) ) )
        }

        # export peak lists

        message( "Info: Exporting the binding site list..." )

        switch( type,
            "gff" = {
                .exportGFF( chrVec=chrVec, muVec=muVec,
                    strengthVec=strengthVec, psize=psize,
                    nameVec=nameVec, filename=filename )
            },
            "bed" = {
                .exportBED( chrVec=chrVec, muVec=muVec,
                    strengthVec=strengthVec, psize=psize,
                    nameVec=nameVec, filename=filename )
            },
            "txt" = {
                .exportTXT( chrVec=chrVec, muVec=muVec,
                    strengthVec=strengthVec, psize=psize,
                    nameVec=nameVec, filename=filename )
            }
        )
    }
)
