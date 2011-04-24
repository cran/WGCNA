#read in prediction fx. It was named version C, now renamed because we'll only use this fx later.

AFcorMI=function(r,m) {
	D=0.43*m^(-0.3)
	epsilon=D^2.2
	out=log(1+epsilon-r^2)/log(epsilon)*(1-D)+D
	out
}




