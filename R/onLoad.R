# Define show methods for large vectors and data.frames
.onLoad <- function(libname, pkgname){
    methods::setClass("advChar",
         slots = list(
           values = "character"
         ))
    methods::setClass("advList",
         slots = list(
           values = "list"
         ))
    setMethod("show", "advChar",
          function(x) {
              if ( length(x@values)<10){
                  print(x@values)
              } else {
                  print(paste0("Vector with ", length(x@values), " values:"))
                  if (is.null(names(x@values))){
                      print(paste0(paste(x@values[1:10], collapse=", "), " ..."))
                  } else {
                      print(paste0(paste(paste0(names(x@values[1:5]),": ",x@values[1:5]),
                                         collapse=", "), " ..."))
                  }
              }
          })
    setMethod("show", "advList",
          function(x) {
              if ( length(x@values)<10){
                  print(x@values)
              } else {
                  print(paste0("Vector with ", length(x@values), " values:"))
                  if (is.null(names(x@values))){
                      print(paste0(paste(x@values[1:10], collapse=", "), " ..."))
                  } else {
                      print(paste0(paste(paste0(names(x@values[1:5]),": ",x@values[1:5]),
                                         collapse=", "), " ..."))
                  }
              }
          })
}