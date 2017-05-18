# S3 classes carry details of output files
XFileOut <- function(format, out.name, out.path = "") {
  if (!is.character(format) || trimws(format) == "") {
    stop("Invalid 'out.format'", call. = F)
  }
  
  if (!is.character(out.name) || trimws(out.name) == "") {
    stop("Invalid 'out.name'", call. = F)
  }
  
  if (!is.character(out.path)) {
    stop("Invalid data type of 'out.path'. Must be string", call. = F)
  }
  
  out.path <- trimws(out.path)
  if (out.path != "") {
    if (!dir.exists(out.path)) {
      tryCatch({
        dir.create(out.path, recursive = T)
      }, warning = function(w) {
        stop("Invalid 'out.path'", call. = F)
      }, error = function(e) {
        stop("Invalid 'out.path'", call. = F)
      })
    }
    
    if (!endsWith(out.path, "/") && !endsWith(out.path, "\\")) {
      out.path <- paste0(out.path, "/")
    }
  }
  
  me <- list(
    suffix = trimws(format),
    name = trimws(out.name),
    path = out.path
  )
  
  me$fullname <- paste0(me$path, me$name, ".", me$suffix)
    
  class(me) <- append(class(me), "XFileOut")
  return(me)
}

RDataFileOut <- function(out.name, out.path = "") {
  me <- XFileOut("RData", out.name, out.path)
  
  class(me) <- append(class(me), "RDataFileOut")
  return(me)
}

TsvFileOut <- function(out.name, out.path = "") {
  me <- XFileOut("tsv", out.name, out.path)
  
  class(me) <- append(class(me), "TsvFileOut")
  return(me)
}

PdfFileOut <- function(out.name, out.path = "") {
  me <- XFileOut("pdf", out.name, out.path)
  
  class(me) <- append(class(me), "PdfFileOut")
  return(me)
}

JpgFileOut <- function(out.name, out.path = "") {
  me <- XFileOut("jpg", out.name, out.path)
  
  class(me) <- append(class(me), "JpgFileOut")
  return(me)
}
