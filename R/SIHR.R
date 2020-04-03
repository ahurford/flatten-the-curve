SIHR <- function(t, y, p) {
	with(as.list(c(y, p)), {
		# The x subscript indicates "no social distancing".
		dSx <- -a * c * Sx * Ix
		dIx <- a * c * Sx * Ix - gamma * Ix - sigma * Ix - v * Ix

		# Hospitalized without a capacity cap for reference
		dHx1 <- sigma * Ix - vH * Hx1 - rho * Hx1

		# Cx: Cumulative cases
		dCx <- a * c * Sx * Ix

		# Hcumx: Cumulative hospital admissions
		dHcumx <- sigma * Ix

		# Much better to code the below as min/max. Poor computation times
		# with if else.
		# Hospitalized with a capacity cap
		# with cap on hospital admissions
		dHx <- sigma * min(Ix, H2) - vH * H0 - rho * H0

		# cumulative unmet need
		dUx <- sigma * max(0, (Ix - H2))

		# with social distancing
		dS <- -a * (1 - m2) * c * S * I
		dI <-  a * (1 - m2) * c * S * I - gamma * I - v * I - sigma * I
		dHcum <- sigma * I

		# without a cap on hospital resources
		dH1 <- sigma * I - vH * H1 - rho * H1

		# culmulative cases
		dC <- a * (1 - m2) * c * S * I

		# with cap on hospital admissions
		dH0 <- sigma * min(I, H2) - vH * H0 - rho * H0

		# cumulative unmet need
		dU <- sigma * max(0, (I - H2))

		list(
			c(
				Sx = dSx,
				Ix = dIx,
				Hx1 = dHx1,
				Hx = dHx,
				Ux = dUx,
				Cx = dCx,
				S = dS,
				I = dI,
				H1 = dH1,
				H0 = dH0,
				U = dU,
				C = dC,
				Hcum = dHcum,
				Hcumx = dHcumx
			)
		)
	})
}
