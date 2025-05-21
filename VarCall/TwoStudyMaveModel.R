model
{
    for (m in 1:M) {
        f.bv[m] ~ dt(mu[m],tau[m],tdf) ## model for replicate values
        tau[m] <- prec.reps[study[m]]*pow(gamma[batchM[m]],-2)
        mu[m] <- beta[batchM[m]] + eta[variantM[m]]*gamma[batchM[m]]
    }

    ## batch 1:48 is study 1 (Mayo), batch (49:150) is study 2 (NCI)
    for (b in 1:48) {
        ## batch-specific random effects
        beta[b] ~ dnorm(mu.beta[1], prec.beta[1])
        gamma[b] ~ dlnorm(mu.gamma[1], prec.gamma[1])
    }
    for (b in 49:B) {
        ## batch-specific random effects
        beta[b] ~ dnorm(mu.beta[2], prec.beta[2])
        gamma[b] ~ dlnorm(mu.gamma[2], prec.gamma[2])
    }

    for (v in 1:V) {
        ## variant-specific function effects:
        eta[v] ~ dnorm(alpha[del[v]], prec.eta[1])  ##LDA (see prec updates below)
    }

    for (v in 1:V){
        ##  Prior on Pathogenicity Status, Pr(D)
        pi.del[v] ~ dbeta(a.pp,b.pp)
        P[v,1] <- (1.0 - pi.del[v])
        P[v,2] <- pi.del[v]
        del[v] ~ dcat(P[v,1:2])
    }

    for (k in 1:2){
        ## measurement variability of replicates
        prec.reps[k] <- pow(sigma2.reps[k], -1.0)
        sigma2.reps[k] <- pow(sigma.reps[k], 2.0)
        sigma.reps[k] ~ dnorm(0.0, cauch.prec.reps[k])
        cauch.prec.reps[k] ~ dgamma(0.5, 0.5)
        ## variability of batch coefficients in mean regression
        prec.beta[k] <- pow(sigma2.beta[k], -1.0)
        sigma2.beta[k] <- pow(sigma.beta[k], 2.0)
        sigma.beta[k] ~ dnorm(0.0, cauch.prec.b[k])
        cauch.prec.b[k] ~ dgamma(0.5, 0.5)
        ## variability of batch scaling in mean regression
        prec.gamma[k] <- pow(sigma2.gamma[k], -1.0)
        sigma2.gamma[k] <- pow(sigma.gamma[k], 2.0)
        sigma.gamma[k] ~ dnorm(0.0, cauch.prec.g[k])
        cauch.prec.g[k] ~ dgamma(0.5, 0.5)
        ## Prior on batch location and scale random effects means
        mu.beta[k] ~ dunif(-10,10)
        mu.gamma[k] ~ dnorm(0.0,1.0)
    }
    ## variability of variant coefficients in mean regression
    for (k in 1:2){
        prec.eta[k] <- pow(sigma2.eta[k], -1.0)
        sigma2.eta[k] <- pow(sigma.eta[k], 2.0)
    }
 
    ## function component variances
    sigma.eta[1] ~ dnorm(0.0, prec.se1)
    sigma.eta[2] <- sigma.eta[1] ##equal
    prec.se1 ~ dgamma(0.5, 0.5)
    prec.se2 ~ dgamma(0.5, 0.5)
    
    prob <- mean(del[])-1
        
}

