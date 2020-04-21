# GOF plot

.plotBIC <- function( object ) {

    maxComp <- get_maxComp(object)

    # generate matrix of BIC & AIC

    bicMat <- aicMat <- matrix( NA, length(get_fits(object)), maxComp )

    for ( i in seq_len(length(get_fits(object))) ) {
        if ( any(is.na(get_optMu(object)[[i]])) ) {
            # error treatment: skip peaks with no fragments
        } else {
            bicMat[i,] <- get_bicVec(object)[[i]]
            aicMat[i,] <- get_aicVec(object)[[i]]
        }
    }

    # GOF plot

    for ( i in seq_len(length(get_stackedFragment(object))) ) {
        # error treatment: skip peaks with no fragments

        # flag if there are no reads in the peak region

        if ( is.na(get_fragSet(object)[[i]][1,1]) ) {
            par( mfrow=c(1,1) )

            plot_title <- paste(get_peakChr(object)[i],": ",
                get_peakStart(object)[i],"-",get_peakEnd(object)[i],sep="")

            plot( 0, 0, type="n", xlab="", ylab="", axes=FALSE,
                main=plot_title, xlim=c(-5,5), ylim=c(-5,5) )
            text( 0, 0, "No fragment (read) in the peak region" )

            next;
        }

        # plot

        par( mfrow=c(1,2) )

        plot_title <- paste("BIC\n",
            get_peakChr(object)[i],": ",get_peakStart(object)[i],"-",get_peakEnd(object)[i],sep="")
        plot( seq_len(maxComp), bicMat[i,],
            xlab="Number of components", ylab="BIC value", main=plot_title )
        xsub <- c(seq_len(maxComp))[ !is.na(bicMat[i,]) ]
        ysub <- bicMat[i,][ !is.na(bicMat[i,]) ]
        lines( xsub, ysub, type="l", lty=2 )

        plot_title <- paste("AIC\n",
            get_peakChr(object)[i],": ",get_peakStart(object)[i],"-",get_peakEnd(object)[i],sep="")
        plot( seq_len(maxComp), aicMat[i,],
            xlab="Number of components", ylab="AIC value", main=plot_title )
        xsub <- c(seq_len(maxComp))[ !is.na(aicMat[i,]) ]
        ysub <- aicMat[i,][ !is.na(aicMat[i,]) ]
        lines( xsub, ysub, type="l", lty=2 )
    }
}
