MapAlignmentAlgorithmIdentification$set("public","align",function(ids,trafos,ref_index)={
  if(!( is_list(ids) && all(sapply(ids,function(i) is_list(i) && all(sapply(i,function(i1) is.R6(i1) && class(i1)[1]=="PeptideIdentification")))) )) { stop("arg id wrong type") }
  if( !(is_list(trafos) && all(sapply(trafos,function(t) is.R6(t) && class(t)[1]=="TransformationDescription"))) ) { stop("arg trafos wrong type") }
  if(!(length(ref_index)==1 && isTRUE(all.equal(ref_index,as.integer(ref_index))))) { stop("arg ref_index wrong type") }
  ids_1 <- r_to_py(ids)
  trafos_1 <- r_to_py(trafos)
  private$py_obj$align(ids_1,trafos_1,as.integer(ref_index))
  ids_1 <- modify_depth(py_to_r(ids_1), 2, function(i) PeptideIdentification$new(i))
  trafos_1 <- modify_depth(py_to_r(trafos_1),1, function(i) TransformationDescription$new(i))
  tryCatch({
          eval.parent(substitute(trafos <- trafos_1))
          eval.parent(substitute(ids <- ids_1))
          invisible()
           }, error = function(e) { invisible() })
} )