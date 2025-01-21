*! version 1.0.0  21jan2025 pdavidboll
version 14
mata

real scalar get_ha_parm_I0(real matrix om_ho, real matrix om_i0, real matrix om_bm, real matrix e) {
    // Variable declarations
    real scalar g, pow, gtry, g1, g2
    real scalar ii

    // Initial search for lower bound where power < 0.5
    pow = 1
    gtry = 1
    while (pow > 0.5) {
        g = gtry;
        // Verify that power is less than 50%
        pow = getpow_qf(om_ho, om_i0 + g * om_bm, e)
        gtry = g / 2
    }
    g1 = g;

    // Initial search for upper bound where power > 0.5
    pow = 0
    gtry = 30
    while (pow < 0.5) {
        g = gtry
        // Verify that power is greater than 50%
        pow = getpow_qf(om_ho, om_i0 + g * om_bm, e)
        gtry = g * 2
    }
    g2 = g

    // Binary search for precise value where power is approximately 0.5
    ii = 1
    while (abs(pow - 0.5) > 0.01) {
        g = (g1 + g2) / 2
        pow = getpow_qf(om_ho, om_i0 + g * om_bm, e)
        if (pow > 0.5) {
            g2 = g
        } else if (pow < 0.5) {
            g1 = g
        }
        ii = ii + 1
        if (ii > 20) {
            break
        }
    }

    return(g)
}

end