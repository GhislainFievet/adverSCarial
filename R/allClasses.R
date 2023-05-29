# Define show methods for large vectors and data.frames

methods::setClass("advChar",
        slots = list(
        values = "character"
        ))
methods::setClass("advList",
        slots = list(
        values = "list"
        ))
setMethod("show", "advChar",
        function(object) {
            if ( length(object@values)<10){
                print(object@values)
            } else {
                print(paste0("Vector with ", length(object@values), " values:"))
                if (is.null(names(object@values))){
                    print(paste0(paste(object@values[1:10], collapse=", "), " ..."))
                } else {
                    print(paste0(paste(paste0(names(object@values[1:5]),": ",object@values[1:5]),
                                        collapse=", "), " ..."))
                }
            }
        })
setMethod("show", "advList",
        function(object) {
            if ( length(object@values)<10){
                print(object@values)
            } else {
                print(paste0("Vector with ", length(object@values), " values:"))
                if (is.null(names(object@values))){
                    print(paste0(paste(object@values[1:10], collapse=", "), " ..."))
                } else {
                    print(paste0(paste(paste0(names(object@values[1:5]),": ",object@values[1:5]),
                                        collapse=", "), " ..."))
                }
            }
        })
