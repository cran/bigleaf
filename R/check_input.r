#######################################
### helper functions to check input ###
#######################################

#' Check Input for Functions in the bigleaf Package
#' 
#' @description Checks length and type of the provided input variables.
#' 
#' @param data   Input data.frame or matrix
#' @param ...    Input variables. Either a list or individual vectors
#' 
#' @note This function can be run for named variables (in which case the return
#'       value will be named according to the given name), or for placeholder
#'       variables that are assigned and named according to e.g. entries of a character 
#'       vector. In the latter case, the input variable has to be named as \code{"var"} or
#'       \code{"var_qc"}.
#' 
#' @keywords internal
check.input <- function(data,...){

  vars <- check.length(list(...))
  
  if (missing(data)){
    data <- NULL
  }
  
  varlist  <- match.call()[-c(1:2)]
  varnames <- c(unlist(sapply(varlist,as.character)))
  varnames <- varnames[!varnames %in% c("c","list")]

  for (v in seq_along(vars)){
    
    var     <- vars[[v]]
    varname <- ifelse(varnames[v] %in% c("var","var_qc"),gsub("\"","",deparse(substitute(var))),varnames[v])
   
    if (is.character(var)){
      if (!missing(data) & !is.null(data)){
        if (length(var) == 1){
          if (var %in% colnames(data)){
            var <- data[,var]
            if (is.numeric(var)){
              assign(varname,var,pos=sys.frame(-1))
            } else {
              stop("column representing '",varname,"' in the input matrix/data.frame must be numeric",call.=FALSE)
            }
          } else {
            stop ("there is no column named '",var,"' in the input matrix/data.frame. Indicate the name of the column representing variable '",varname,"', or alternatively, provide a numeric vector of the same length as the input matrix/data.frame or of length 1.",call.=FALSE)
          }
        } else {
          stop("name of variable '",varname,"' must have length 1",call.=FALSE)
        }
      } else {
        if ("data" %in% names(formals(sys.function(which=-1)))){
          if (var %in% as.character(unlist(match.call(definition=sys.function(-1),call=sys.call(-1))[-1]))){
            stop("variable '",var,"' is of type character and interpreted as a column name, but no input matrix/data.frame is provided. Provide '",var,"' as a numeric vector, or an input matrix/data.frame with a column named '",var,"'",call.=FALSE)
          } else {
            stop("variable '",var,"' is not provided",call.=FALSE)
          }
        } else {
          stop("variable '",var,"' must be numeric",call.=FALSE)
        }
      }
    } else {
      if (length(var) < 2){
        if (is.null(var)){
          assign(varname,var,pos=sys.frame(-1))
          next
        } else if (is.na(var)){
          assign(varname,var,pos=sys.frame(-1))
          next
        }
      }
      if (!missing(data) & !is.null(data)){
        if (is.numeric(var) & length(var) == nrow(data)){
          assign(varname,var,envir=sys.frame(-1))
        } else if (is.numeric(var) & length(var) != nrow(data)) {
          if (length(var) == 1){
            var <- rep(var,length=nrow(data))
            assign(varname,var,envir=sys.frame(-1))
          } else {
            stop("variable '",varname,"' must have the same length as the input matrix/data.frame or length 1. Do NOT provide an input matrix/data.frame if none of its variables are used!",call.=FALSE)
          }
        } else if (!is.numeric(var)){
          stop("variable '",varname,"' must be numeric",call.=FALSE)
        }
      } else {
        if (is.numeric(var)){
          assign(varname,var,envir=sys.frame(-1))
        } else {
          stop("variable '",varname,"' must be numeric",call.=FALSE)
        }
      }
    }
  }
}



#' Test Variables for Equal Length
#' 
#' @param varlist List of variables for which the length has to be compared
#' 
#' @note This function only plays a role if no input data.frame or matrix are 
#'       provided. In this case it ensures that provided vectors have the same
#'       length to avoid trouble later up the function call.
#'       
#' @keywords internal
check.length <- function(varlist){
  
  if (is.list(unlist(varlist,recursive=FALSE))){
    varlist <- unlist(varlist,recursive=FALSE)
  }
  
  length.vars <- sapply(varlist,length)
  length.vars <- length.vars[length.vars > 0]
  
  if (length(unique(length.vars)) >= 2){
    if (sort(unique(length.vars))[1] != 1 | length(unique(length.vars)) > 2){
      stop("All input variables must have the same length or a length of 1!",call.=FALSE)
    }
  }
  return(varlist)
}
