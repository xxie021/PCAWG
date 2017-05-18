source("source/transform.R")
source("source/utils.R")

# Saves a motif matrix as a RData or TSV file with the specific format. If
# {@code file.out} is not specified, prints the result to the STDOUT.
# @param  mut.ctx         a motif matrix
# @param  type            one of the five types accepted by the
#                         {@code Transform3} function:
#                         (1) base.summary
#                         (2) base.spectrum
#                         (3) base.ctx.summary
#                         (4) base.ctx.heatmap.plot
#                         (5) base.ctx.heatmap.real
# @param  file.out        a subtype of the {@code XFileOut} object (currently
#                         accepts {@code RDataFileOut} or {@code TsvFileOut}
#                         only)
# @param  rdata.obj.name  the customised variable name of an RData object
Summary <- function(mut.ctx, type, file.out = NULL, rdata.obj.name = NULL) {
  data <- Transform3(mut.ctx, type)
  
  if (is.object(file.out) && class(file.out)[2] == "XFileOut") {
    if (class(file.out)[3] == "RDataFileOut") {
      if (is.character(rdata.obj.name)) {
        Save2RData(data, file.out, rdata.obj.name)
      } else {
        Save2RData(data, file.out)
      }
    } else if (class(file.out)[3] == "TsvFileOut") {
      cat("Info: Saving file \"", basename(file.out$fullname), "\" ...\n",
          sep = "")
      write.table(data, file.out$fullname, quote = F, sep = "\t",
                  col.names = NA)
      cat("Info: File \"", basename(file.out$fullname), "\" saved\n", sep = "")
    } else {
      stop("Output format is not currently supported", call. = F)
    }
  } else {
    print(data)
  }
}
