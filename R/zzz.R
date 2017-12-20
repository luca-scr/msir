.onAttach <- function(lib, pkg)
{
  version <- packageVersion("msir")
  packageStartupMessage("Package 'msir' version ", version)
  packageStartupMessage("Type 'citation(\"msir\")' for citing this R package in publications.")
  invisible()
}
  