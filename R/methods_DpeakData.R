
# generic methods for "DpeakData" class

setMethod(
  f = "get_fragSet",
  signature = "DpeakData",
  definition = function(x) x@fragSet
)

setMethod(
  f = "get_PET",
  signature = "DpeakData",
  definition = function(x) x@PET
)

setMethod(
  f = "get_fragLenTable",
  signature = "DpeakData",
  definition = function(x) x@fragLenTable
)

setMethod(
  f = "get_Fratio",
  signature = "DpeakData",
  definition = function(x) x@Fratio
)

setMethod(
  f = "get_aveFragLen",
  signature = "DpeakData",
  definition = function(x) x@aveFragLen
)

setMethod(
  f = "get_stackedFragment",
  signature = "DpeakData",
  definition = function(x) x@stackedFragment
)

setMethod(
  f = "get_peakChr",
  signature = "DpeakData",
  definition = function(x) x@peakChr
)

setMethod(
  f = "get_peakStart",
  signature = "DpeakData",
  definition = function(x) x@peakStart
)

setMethod(
  f = "get_peakEnd",
  signature = "DpeakData",
  definition = function(x) x@peakEnd
)

setMethod(
  f = "get_emptyList",
  signature = "DpeakData",
  definition = function(x) x@emptyList
)

setMethod(
    f="show",
    signature="DpeakData",
    definition=function( object ) {
        # extract objects

        peakChr <- get_peakChr(object)
        chrList <- unique(peakChr)
        PET <- get_PET(object)
        aveFragLen <- get_aveFragLen(object)

        # summary

        nFrag <- unlist( lapply( get_fragSet(object),
            function(x) ifelse( !is.null(x), nrow(x), 0 )
         ) )
        sumRead <- sum( nFrag )
        medNumRead <- median( nFrag )

        cat( "------------------------------------------------------------\n" )
        cat( "Summary: Dpeak data (class: DpeakData)\n" )
        cat( "------------------------------------------------------------\n" )
        cat( "Number of peaks: ",length(peakChr),"\n", sep="" )
        cat( "Number of chromosomes: ",length(chrList),"\n", sep="" )
        cat( "Tag type: ",ifelse(PET,"PET","SET"),"\n", sep="" )
        if ( PET ) {
            cat( "Median fragment length: ",aveFragLen,"\n", sep="" )
        } else {
            cat( "Fragment length (provided by user): ",aveFragLen,"\n", sep="" )
        }
        cat( "Number of utilized reads: ",sumRead,"\n", sep="" )
        cat( "Median number of reads in each peak: ",medNumRead,"\n", sep="" )
        cat( "------------------------------------------------------------\n" )
    }
)

setMethod(
    f="print",
    signature="DpeakData",
    definition=function( x ) {
        warning( "'print' method for 'DpeakData' class is not supported yet." )
    }
)

setMethod(
    f="printEmpty",
    signature="DpeakData",
    definition=function( object, ... ) {
        if ( identical(get_emptyList(object)[1],"") ) {
            message( "Every peak region has at least one read." )
        } else {
            out <- apply( as.matrix(get_emptyList(object)), 1,
                function(x) {
                    splitOrg <- strsplit( x, "_" )[[1]]

                    if ( length(splitOrg) > 3 ) {
                        xSplit <- rep( NA, 3 )
                        xSplit[1] <- paste( splitOrg[ seq_len(length(splitOrg)-2) ],
                            collapse="_" )
                        xSplit[seq(from = 2, to = 3)] <- splitOrg[ -c(seq_len(length(splitOrg)-2)) ]
                    } else {
                        xSplit <- splitOrg
                    }
                    return(xSplit)
            } )

            out <- data.frame( t(out), stringsAsFactors=FALSE )
            out[,2] <- as.numeric(out[,2])
            out[,3] <- as.numeric(out[,3])
            colnames(out) <- c( "chrID", "peakStart", "peakEnd" )

            return( out )
        }
    }
)

setMethod(
    f="plot",
    signature=c("DpeakData","missing"),
    definition=function( x, y, filename=NULL, strand=FALSE, extension=1, smoothing=FALSE, ... ) {
        # extract estimates

        PET <- get_PET(x)

        # error treatment

        if ( extension < 1 ) {
            stop( "Error: A negative 'extension' value is not allowed." )
        }

        # plot

        pdf( filename )

        for ( i in seq_len(length(get_stackedFragment(x)) )) {

            plot_title <- paste(get_peakChr(x)[i],": ",
                                get_peakStart(x)[i],"-",
                                get_peakEnd(x)[i],sep="")

            # flag if there are no reads in the peak region

            if ( is.na(get_fragSet(x)[[i]][1,1]) ) {

                plot( 0, 0, type="n", xlab="", ylab="", axes=FALSE,
                    main=plot_title, xlim=c(-5,5), ylim=c(-5,5) )
                text( 0, 0, "No fragment (read) in the peak region" )

                next;
            }

            # plot

            xlim <- rep( NA, 2 )
            xlim[1] <- min( get_peakStart(x)[i], get_stackedFragment(x)[[i]][,1] )
            xlim[2] <- max( get_peakEnd(x)[i], get_stackedFragment(x)[[i]][,1] )

            if ( strand ) {

                .plotStrandData( stackedFragment=get_stackedFragment(x)[[i]],
                    fragSet=get_fragSet(x)[[i]], plot_title=plot_title, xlim=xlim,
                    PET=PET, extension=extension, smoothing=smoothing )

                if ( extension == 1 ) {
                    legend( "topright", lty=c(1,1), col=c("black","red"),
                        c("Stacked reads (forward)","Stacked reads (reverse)")
                    )
                } else {
                    legend( "topright", lty=c(1,1), col=c("black","red"),
                        c("Stacked extended reads (forward)","Stacked extended reads (reverse)")
                    )
                }
            } else {
                # combined

                plot( get_stackedFragment(x)[[i]][,1],
                      get_stackedFragment(x)[[i]][,2], type="l",
                    xlab="Genomic coordinates", ylab="Frequency",
                    main=plot_title,
                    xlim=xlim, ylim=c(0,max(get_stackedFragment(x)[[i]][,2])*1.2) )

                legend( "topright", lty=1, col="black", "Stacked fragments" )
            }
            abline( v=get_peakStart(x)[i], col="red", lty=2 )
            abline( v=get_peakEnd(x)[i], col="red", lty=2 )
        }

        dev.off()
    }
)
