# GOF plot

.plotGOF <- function( object, nsimul=10000, seed=12345,
    extension=1, smoothing=TRUE, nCore=8 ) {

    # extract estimates

    maxComp <- get_maxComp(object)
    peakStart <- get_peakStart(object)
    peakEnd <- get_peakEnd(object)

    PET <- get_PET(object)
    fragLenTable <- get_fragLenTable(object)
    aveFragLen <- get_aveFragLen(object)
    Fratio <- get_Fratio(bject)

    # generate list to generate fragments

    optList <- vector( "list", length(get_optMu(object)) )
    for ( i in seq_len(length(get_optMu(object))) ) {
        # error treatment: skip peaks with no fragments

        if ( is.na(get_fragSet(object)[[i]][1,1]) ) {
            optList[[i]] <- matrix( NA )
        } else {
            # stack info

            optList_i <- list()

            optList_i$nsimul <- nsimul
            optList_i$mu <- get_optMu(object)[[i]]
            optList_i$pi <- get_optPi(object)[[i]]
            optList_i$pi0 <- get_optPi0(object)[[i]]
            optList_i$minS <- min(get_stackedFragment(object)[[i]][,1])
            optList_i$maxS <- max(get_stackedFragment(object)[[i]][,1])
            optList_i$peakstart <- peakStart[i]
            optList_i$peakend <- peakEnd[i]
            if ( PET == FALSE ) {
                optList_i$delta <- get_optDelta(object)[[i]]
                optList_i$sigma <- get_optSigma(object)[[i]]
            }
            optList[[i]] <- optList_i
        }
    }

    if ( PET ) {
        Lvalue <- as.numeric(as.vector(names(fragLenTable)))
        Lprob <- as.numeric(fragLenTable)
    }

    # generate fragments (using parallel computing, if parallel exists)

    message( "Info: Generating simulated fragments from the fitted model..." )

    if ( is.element( "parallel", installed.packages()[,1] ) ) {
        # if "parallel" package exists, utilize parallel computing with "parallel::mclapply"

        simList <- parallel::mclapply( optList, function(x) {
            simFrag <- .generateFragment( object=x,
                PET=PET, Lvalue=Lvalue, Lprob=Lprob,
                Fratio=Fratio, aveFragLen=aveFragLen )
            simStack <- .stackFragment( simFrag )
            out <- list( stackedFragment=simStack, fragSet=simFrag )
            return( out )
        }, mc.cores = nCore )
    } else {
        # otherwise, use usual "lapply"

        simList <- lapply( optList, function(x) {
            simFrag <- .generateFragment( object=x,
                PET=PET, Lvalue=Lvalue, Lprob=Lprob,
                Fratio=Fratio, aveFragLen=aveFragLen )
            simStack <- .stackFragment( simFrag )
            out <- list( stackedFragment=simStack, fragSet=simFrag )
            return( out )
        } )
    }

    # GOF plot

    message( "Info: Generating GOF plots..." )

    for ( i in seq_len(length(get_stackedFragment(object))) ) {

        plot_title <- paste(get_peakChr(object)[i],": ",
            get_peakStart(object)[i],"-",get_peakEnd(object)[i],sep="")

        # flag if there are no reads in the peak region

        if ( is.na(get_fragSet(object)[[i]][1,1]) ) {

            plot( 0, 0, type="n", xlab="", ylab="", axes=FALSE,
                main=plot_title, xlim=c(-5,5), ylim=c(-5,5) )
            text( 0, 0, "No fragment (read) in the peak region" )

            next;
        }

        # plot

        stackedSimFrag <- simList[[i]]$stackedFragment

        xlim <- rep( NA, 2 )
        xlim[1] <- min( get_peakStart(object)[i], get_stackedFragment(object)[[i]][,1],
            stackedSimFrag[,1] )
        xlim[2] <- max( get_peakEnd(object)[i], get_stackedFragment(object)[[i]][,1],
            stackedSimFrag[,1] )

        if ( PET == FALSE ) {
            .plotStrandData( stackedFragment=get_stackedFragment(object)[[i]],
                fragSet=get_fragSet(object)[[i]], plot_title=plot_title, xlim=xlim,
                PET=PET, extension=extension, smoothing=smoothing )
        } else {
            # PET

            plot( get_stackedFragment(object)[[i]][,1],
                get_stackedFragment(object)[[i]][,2], type="l",
                xlab="Genomic coordinates", ylab="Frequency",
                main=plot_title,
                xlim=xlim, ylim=c(0,max(get_stackedFragment(object)[[i]][,2])*1.2) )
        }
        abline( v=get_peakStart(object)[i], col="red", lty=2 )
        abline( v=get_peakEnd(object)[i], col="red", lty=2 )

        # plot simulated data

        if ( PET == TRUE ) {
            # PET

            simX <- stackedSimFrag[,1]
            simY <- stackedSimFrag[,2] * max(get_stackedFragment(object)[[i]][,2]) /
                max(stackedSimFrag[,2])
                # adjust frequency to make frequency comparable to observed one
            lines( simX, simY, col="gray", lty=2 )

            legend( "topright", lty=c(1,2), col=c("black","gray"),
                c("Stacked fragments (data)","Stacked fragments (simulated)") )
        } else {
            # SET

            ratio <- max(get_stackedFragment(object)[[i]][,2]) /
                max(stackedSimFrag[,2])

            .lineStrandData( stackedFragment=stackedSimFrag,
                fragSet=simList[[i]]$fragSet, adjustment=ratio,
                plot_title=plot_title, xlim=xlim,
                PET=PET, extension=extension, smoothing=smoothing )

            if ( extension == 1 ) {
                legend( "topright", lty=c(1,1,2,2),
                    col=c("black","red","gray","pink"),
                    c("Stacked reads: forward (data)",
                    "Stacked reads: reverse (data)",
                    "Stacked reads: forward (simulated)",
                    "Stacked reads: reverse (simulated)")
                )
            } else {
                legend( "topright", lty=c(1,1,2,2),
                    col=c("black","red","gray","pink"),
                    c("Stacked extended reads: forward (data)",
                    "Stacked extended reads: reverse (data)",
                    "Stacked extended reads: forward (simulated)",
                    "Stacked extended reads: reverse (simulated)")
                )
            }
        }
    }

    message( "Info: Done!" )
}
