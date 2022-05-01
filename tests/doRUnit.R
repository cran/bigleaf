#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++ call testthat unit tests in directory inst\tests +++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

if (require("testthat", quietly = TRUE)) {
	pkg <- "bigleaf"
	#library(pkg, character.only = TRUE)
	#test_package(pkg)
	test_check(pkg)
} else {
	warning("cannot run unit tests -- package testthat is not available")
}

