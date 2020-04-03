SIR <- function(t, y, p) {
	with(as.list(c(y, p)), {
		# No social distancing
		# -----------------
		# Sx: Susceptible (no social distancing)
		dSx <- -a * c * Sx * Ix
		# Ix: Infected (no social distancing)
		dIx <- a * c * Sx * Ix - gamma * Ix - v * Ix
		# Fx: Cumultative fatalities (no)
		dFx <- v * Ix
		# Cumulative cases
		dCx <- a*c*Sx*Ix
		# Variables as above but with social distancing
		# Social distancing is m1
		dS <- -a * (1 - m1) * c * S * I
		dI <-  a * (1 - m1) * c * S * I - gamma * I - v * I
		dF <- v * I
		dC <- a*c*(1-m1)*S*I

		list(c(Sx = dSx, Ix = dIx, Fx = dFx, Cx = dCx, S = dS,
					 I = dI, Fs = dF, C=dC))
	})
}
