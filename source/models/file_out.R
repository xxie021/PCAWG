# S3 classes carry details of output files
# Four inherited classes: "RDataFileOut", "TsvFileOut", "PdfFileOut" and
# "JpgFileOut"
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
      dir.create(out.path)
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

GenFullFileName <- function(file.out) {
  UseMethod("GenFullFileName", file.out)
}

GenFullFileName.XFileOut <- function(file.out) {
  return(paste0(file.out$path, file.out$name, ".", file.out$suffix))
}
