
# generic methods for "DpeakMotif" class

setMethod(
  f = "get_motif",
  signature = "DpeakMotif",
  definition = function(x) x@motif
)

setMethod(
  f = "get_locMotif",
  signature = "DpeakMotif",
  definition = function(x) x@locMotif
)

setMethod(
    f="show",
    signature="DpeakMotif",
    definition=function( object ) {
      # extract objects

      motifIden <- get_motif(object)
      locMotif <- get_locMotif(object)

      if ( all(is.na(unlist(locMotif))) ) {
      	# if no peak includes a motif

      	cat( "------------------------------------------------------------\n" )
      	cat( "Info: Preprocessing summary\n" )
      	cat( "------------------------------------------------------------\n" )
      	cat( "Identified motif:\n", sep="" )
      	for ( i in seq_len(length(motifIden)) ) {
      		cat( "\t",motifIden[i],"\n", sep="" )
      	}
      	cat( "Note: No peak contains detected motifs.\n")
      	cat( "------------------------------------------------------------\n" )
      } else {
		    # info about preprocessing
		    # - identified motif, # detected motifs

      	listMotifVec <- list()
      	k <- 1
      	for ( i in seq_len(length(locMotif)) ) {
      		if ( !is.na(locMotif[[i]][1]) ) {
      			listMotifVec[[k]] <- locMotif[[i]]
      			k <- k + 1
      		}
      	}

		    cat( "------------------------------------------------------------\n" )
		    cat( "Info: Preprocessing summary\n" )
		    cat( "------------------------------------------------------------\n" )
		    cat( "Identified motif:\n", sep="" )
		    for ( i in seq_len(length(motifIden)) ) {
		    	cat( "\t",motifIden[i],"\n", sep="" )
		    }
		    cat( "Number of peaks containing detected motifs: ",length(listMotifVec),"\n", sep="" )
		    cat( "Number of peaks missing detected motifs: ",( length(locMotif) - length(listMotifVec) ),"\n", sep="" )
		    cat( "Median number of regulatory sequences in each peak: ",median( sapply( listMotifVec, length ) ),"\n", sep="" )
		    cat( "\t(after excluding peaks missing detected motifs)\n" )
		    cat( "------------------------------------------------------------\n" )
    	}
    }
)
