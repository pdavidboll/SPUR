*! version 1.0.0  21jan2025 pdavidboll
version 14
mata

real matrix get_sigma_dm(real matrix distmat, real scalar c) {
    // Variable declarations
    real scalar n
    real matrix sigma, sigma_dm

    // Step 1: Compute sigma as exp(-c * distmat)
    n = rows(distmat)
    sigma = exp(-c * distmat)

    // Step 2: Demean sigma
    sigma_dm = demeanmat(sigma)

    return(sigma_dm)
}

end