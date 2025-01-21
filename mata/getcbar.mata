*! version 1.0.0  21jan2025 pdavidboll
version 14
mata

real scalar getcbar(real scalar rhobar, real matrix distmat) {
	real scalar c0
	real scalar c1
	real vector vd
	real scalar v, i1
	real scalar jj
	real scalar cm, cbar
	
	c0 = 10
	c1 = 10
	vd = lvech(distmat)

    // Adjust c0 and c1 bounds so c0 yields v > rhobar and c1 yields v < rhobar
    i1 = 0
	jj = 0

    while (i1 == 0) {
        v = mean(exp(-c0 * vd))
        i1 = (v > rhobar)
      
        if (i1 == 0) {
            c1 = c0
            c0 = c0 / 2
            jj = jj + 1
        }
        if (jj > 500) {
            _error("rhobar too large")
        }
    }

    // Verify that c1 yields a value less than rhobar
    i1 = 0
    jj = 0
    while (i1 == 0) {
        v = mean(exp(-c1 * vd))
        i1 = (v < rhobar)
      
        if (i1 == 0) {
            c0 = c1
            c1 = 2 * c1
            jj = jj + 1
        }
        if (c1 > 10000) {
            i1 = 1
        }
        if (jj > 500) {
            _error("rhobar too small")
        }
    }

    // Bisection determination of cbar
    while ((c1 - c0) > 0.001) {
        cm = sqrt(c0 * c1)
        v = mean(exp(-cm * vd))
      
        if (v < rhobar) {
            c1 = cm
        } else {
            c0 = cm
        }
    }

    cbar = sqrt(c0 * c1)
    return(cbar)
}

end