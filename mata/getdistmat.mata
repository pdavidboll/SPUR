*! version 1.0.0  21jan2025 pdavidboll
version 14
mata

real matrix getdistmat(real matrix s)
// computes matrix of distances from locations
// if external latlongflag=0, distances computed as norm, otherwise latitude / longitude and haversine formula for sphere with radius 1/Pi
{	
	external real scalar latlongflag
	real matrix mat
	real scalar n,i,c
	real matrix d
	n=rows(s)
	mat=J(n,n,0)
	if(latlongflag==0){
		for(i=1;i<=n;i++){
			mat[.,i]=sqrt(rowsum((s:-s[i,.]):^2))
		}
	}
	else
		{
		c=3.14159265359/180
		for(i=1;i<=n;i++){
			d=(.5*c)*(s:-s[i,.])
			mat[.,i]=asin(sqrt(sin(d[.,1]):^2  + cos(c*s[i,1])*(cos(c*s[.,1]):*(sin(d[.,2]):^2))))/3.14159265359
		}
	}	
	return(mat)
}	

end