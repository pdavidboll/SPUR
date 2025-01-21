*! version 1.0.0  21jan2025 pdavidboll
version 14
mata

real scalar getpow_qf(real matrix om0, real matrix om1, real matrix e) {
    // Variable declarations
    real matrix om0i, om1i, ch_om0, ch_om1, ch_om0i, ch_om1i
    real matrix ho, ha, ya_o, yo_a
    real vector qe, qa_o, qo_a, lr_o, lr_a
    real scalar cv, pow

    // Step 1: Compute inverses
    om0i = invsym(om0)    // Use invsym() for symmetric matrix inverse
    om1i = invsym(om1)
  
    // Step 2: Cholesky decompositions
    ch_om0 = cholesky(om0)'
    ch_om1 = cholesky(om1)'
    ch_om0i = cholesky(om0i)'
    ch_om1i = cholesky(om1i)'
  
    // Step 3: Compute transformed matrices
    ho = ch_om1i * ch_om0'
    ha = ch_om0i * ch_om1'
  
    // Step 4: Compute quadratic forms
    qe = colsum(e:^2)    // Sum of squares of columns of e
    ya_o = ho * e
    yo_a = ha * e
    qa_o = colsum(ya_o:^2) // Sum of squares of columns of ya_o
    qo_a = colsum(yo_a:^2) // Sum of squares of columns of yo_a
  
    // Step 5: Likelihood ratios
    lr_o = qe :/ qa_o
    lr_a = qo_a :/ qe
  
    // Step 6: Compute 95th percentile critical value and calculate power
    cv = mm_quantile(lr_o', 1, 0.95) // 95th percentile of lr_o
    pow = mean(lr_a' :> cv)
  
    return(pow)
}

end