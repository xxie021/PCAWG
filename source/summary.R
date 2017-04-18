if (!exists("Transform3", mode = "function")) {
  source("source/transform.R")
}

if (!exists("Save2RData", mode = "function")) {
  source("source/utils.R")
}

if (!exists("XFileOut", mode = "function")) {
  source("source/models/file_out.R")
}

# This function calls the function "Transform3", so both "mut.ctx" and "type"
# share the same constraints of that one. Currently it accepts two types of
# "file.out": "RDataFileOut" and/or "TsvFileOut". It also accepts an optional
# parameter "rdata.obj.name" for saving the output file in the "RData" format
# using a different name
Summary <- function(mut.ctx, type, file.out = NA, ...) {
  data <- Transform3(mut.ctx, type)
  
  if (is.object(file.out) && class(file.out)[2] == "XFileOut") {
    if (file.out$suffix == "RData") {
      args <- list(...)
      if ("rdata.obj.name" %in% names(args)) {
        Save2RData(data, file.out, args$rdata.obj.name)
      } else {
        Save2RData(data, file.out)
      }
    } else if (file.out$suffix == "tsv") {
      file.fullname <- GenFullFileName(file.out)
      cat("Info: Saving file \"", basename(file.fullname), "\" ...\n", sep = "")
      write.table(data, file.fullname, quote = F, sep = "\t", col.names = NA)
      cat("Info: File \"", basename(file.fullname), "\" saved\n", sep = "")
    } else {
      stop("Unrecognised output format", call. = F)
    }
  } else {
    print(data)
  }
}
