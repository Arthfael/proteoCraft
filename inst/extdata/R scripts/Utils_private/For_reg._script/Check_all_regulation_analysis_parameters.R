options(stringsAsFactors = FALSE)
setwd("D:\\groups_temp\\")
require(proteoCraft)
require(magrittr)

fls <- list.files(recursive = TRUE)
params <- grep("combined/txt/Parameters\\.csv$|combined/Parameters\\.csv$", fls, value = TRUE)
params2 <- lapply(params, function(x) {
  Param.load(x, filter.deprecated = TRUE)
})
param_names <- unique(unlist(sapply(params2, colnames)))
param_descript <- sapply(param_names, function(x) {
  sapply(1:length(params), function(y) {
    
  })
})
