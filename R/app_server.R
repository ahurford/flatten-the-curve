# App Server

#' @import shiny
app_server <- function(input, output,session) {

	# Data NL
	callModule(mod_data_NL_server, "data_NL_1")

}
