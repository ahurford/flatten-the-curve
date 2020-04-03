SEIR <- function(t, y, p) {
	with(as.list(c(y, p)), {

		if (t <= today) {
			dS <- -beta * S * (I + E)
			dE <- beta * S * (I + E) - a1 * E
			dI <- a1 * E - gamma1 * I - v1 * I
			dC <- a1 * E
		} else {
			dS <- -beta1 * S * (I + E)
			dE <- beta1 * S * (I + E) - a1 * E
			dI <- a1 * E - gamma1 * I - v1 * I
			dC <- a1 * E
		}

		list(c(
			S = dS,
			E = dE,
			I = dI,
			C = dC
		))
	})
}


SEIRnull <- function(t, y, p) {
	with(as.list(c(y, p)), {

		dS <- -beta * S * (I + E)
		dE <- beta * S * (I + E) - a1 * E
		dI <- a1 * E - gamma1 * I - v1 * I
		dC <- a1 * E

		list(c(
			S = dS,
			E = dE,
			I = dI,
			C = dC
		))
	})
}
