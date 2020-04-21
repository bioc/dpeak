# plot fitting results

.plotFit <- function( object, threshold=1000,
    strand=FALSE, extension=1, smoothing=FALSE ) {

    # extract estimates

    PET <- get_PET(object)
    psize <- get_psize(object)

    # error treatment

    if ( extension < 1 ) {
        stop( "Negative 'extension' is not allowed!" )
    }

    # plot

    for ( i in seq_len(length(get_stackedFragment(object))) ) {

        plot_title <- paste(get_peakChr(object)[i],": ",get_peakStart(object)[i],"-",get_peakEnd(object)[i],sep="")

        # flag if there are no reads in the peak region

        if ( is.na(get_fragSet(object)[[i]][1,1]) ) {

            plot( 0, 0, type="n", xlab="", ylab="", axes=FALSE,
                main=plot_title, xlim=c(-5,5), ylim=c(-5,5) )
            text( 0, 0, "No fragment (read) in the peak region" )

            next;
        }

        # plot

        xlim <- rep( NA, 2 )
        xlim[1] <- min( get_peakStart(object)[i], get_stackedFragment(object)[[i]][,1] )
        xlim[2] <- max( get_peakEnd(object)[i], get_stackedFragment(object)[[i]][,1] )

        Ni <- nrow(get_fragSet(object)[[i]])

        if ( strand==TRUE ) {

            .plotStrandData( stackedFragment=get_stackedFragment(object)[[i]],
                fragSet=get_fragSet(object)[[i]], plot_title=plot_title, xlim=xlim,
                PET=PET, extension=extension, smoothing=smoothing )

            if ( extension == 1 ) {
                legend( "topright", lty=c(1,1,2,2), lwd=c(1,1,2,2),
                    col=c("black","red","blue","lightblue"),
                    c("Stacked reads (forward)","Stacked reads (reverse)",
                        paste("Binding sites (strength > ",threshold,")",sep=""),
                        paste("Binding sites (strength <= ",threshold,")",sep=""))
                )
            } else {
                legend( "topright", lty=c(1,1,2,2), lwd=c(1,1,2,2),
                    col=c("black","red","blue","lightblue"),
                    c("Stacked extended reads (forward)","Stacked extended reads (reverse)",
                        paste("Binding sites (strength > ",threshold,")",sep=""),
                        paste("Binding sites (strength <= ",threshold,")",sep=""))
                )
            }
        } else {
            # combined

            plot( get_stackedFragment(object)[[i]][,1], get_stackedFragment(object)[[i]][,2], type="l",
                xlab="Genomic coordinates", ylab="Frequency",
                main=plot_title,
                xlim=xlim, ylim=c(0,max(get_stackedFragment(object)[[i]][,2])*1.2) )

            legend( "topright", lty=c(1,2,2), lwd=c(1,2,2),
                col=c("black","blue","lightblue"),
                c("Stacked fragments",
                    paste("Binding sites (strength > ",threshold,")",sep=""),
                    paste("Binding sites (strength <= ",threshold,")",sep="")
                )
            )
        }

        abline( v=get_optMu(object)[[i]][ Ni * get_optPi(object)[[i]] > threshold ], col="blue", lty=2, lwd=2 )
        abline( v=get_optMu(object)[[i]][ Ni * get_optPi(object)[[i]] <= threshold ], col="lightblue", lty=2, lwd=2 )
        abline( v=get_peakStart(object)[i], col="red", lty=2 )
        abline( v=get_peakEnd(object)[i], col="red", lty=2 )
    }
}
