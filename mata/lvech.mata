*! version 1.0.0  21jan2025 pdavidboll
version 14
mata

real vector lvech(real matrix S) {
    real scalar n
	real vector v
	real scalar idx, i, j
	
	n = rows(S)
	idx = 1

    // Preallocate v with the correct size for efficiency
    v = J((n * (n - 1)) / 2, 1, .)

    // Collect elements below the main diagonal
	for (j = 1; j <= n; j++) {
		for (i = j+1; i <= n; i++) {
            v[idx] = S[i, j]
            idx++
        }
    }
    
    return(v)
}

end