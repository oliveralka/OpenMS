Attachment$set("active","tableRows", function(tableRows) = {
  if(!missing(tableRows)){
    if( !(is_list(tableRows) && lapply(tableRows,function(v0) (is_list(v0) || is_vector(v0)) && sapply(v0,function(v0_1) (is_scalar_character(v0_1) || is.R6(v0_1) && class(v0_1)[1]=="String")))) )
      {
          private$py_obj$tableRows <- tableRows
      }
  } else {
        return(private$py_obj$tableRows)
  }
} )

