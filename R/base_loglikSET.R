
# log likelihood

.loglikSET <- function( S, E, strand, mu, pi, pi0,
    delta, sigma, Fratio, sindex, beta, R, alpha ) {

    # initialization

    FRvec <- sindex$FRvec
    indF <- sindex$indF
    indR <- sindex$indR

    N <- length(S)
    n_group <- length(mu)

    out <- alpha

    # likelihood: log P(L) is not incorporated (not affect loglik)

    out_g <- rep( 0, N )

    for ( g in seq_len(n_group) ) {
        termF <- Fratio * dnorm( S, mu[g] - delta, sigma ) * indF
        termR <- ( 1 - Fratio ) * dnorm( E, mu[g] + delta, sigma ) * indR

        out_g <- out_g + pi[g] * ( termF + termR )
    }
    termBG <- FRvec * pi0 / ( R + beta - 1 )

    out_g <- out_g + termBG

    nonzero <- which( out_g > 0 )
    out <- sum(log( ( out * out_g )[ nonzero ] ))

    return(out)
}
