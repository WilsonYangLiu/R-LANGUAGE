Pvalue.norm.sim <- function(n=50, mu=0, mu0=0, sigma=1, sigma0=sigma,
                            test=c('z','t'),
                            alternative=c('two.sided', 'less', 'greater',
                            '<>','!=','<','>'),
                            alpha=0.05, B=1e4) {
    test <- match.arg(test)
    alternative <- match.arg(alternative)

    x <- matrix(rnorm( n*B, mu, sigma ), nrow=n)
    xbar <- colMeans(x)

    if( is.na(sigma0) ) sigma0 <- apply(x, 2, sd)

    ts <- (xbar - mu0)/sigma0*sqrt(n)

    pdist <- switch(test,
               z=function(x, lower.tail) pnorm(x, lower.tail=lower.tail),
               t=function(x, lower.tail) pt(x, df=n-1, lower.tail=lower.tail)
                    )

    p.vals <- switch(alternative,
                     '!='=,'<>'=,
                     two.sided = 2*pmin( pdist(ts,TRUE), pdist(ts,FALSE) ),
                     '<'=,
                     less      = pdist(ts, TRUE),
                     '>'=,
                     greater   = pdist(ts, FALSE)
                     )

    op <- par(mfrow=c(2,1))
    hist(p.vals, main='', xlab='P-Values')

    if( !is.na(alpha) ) {
        abline(v=alpha, col='red')
        title(sub=paste( round(mean(p.vals <= alpha)*100, 1), '% <= ',
              alpha))
    }

    qqplot( seq(along=p.vals)/(B+1), p.vals,
           xlab='Theoretical quantiles of Uniform',
           ylab='P-values')
    abline(0,1, col='grey')

    par(op)

    invisible(p.vals)
}
