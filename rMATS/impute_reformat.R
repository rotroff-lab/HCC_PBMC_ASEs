#!/usr/bin/env Rscript

library(optparse)
library(parallel)

# Command-line options
option_list <- list(
  make_option(c("-s", "--srp_file"), type="character", help="SRP file list", metavar="character"),
  make_option(c("-m", "--metadata_path"), type="character", help="Path to metadata", metavar="character"),
  make_option(c("-c", "--comparison"), type="character", help="Comparison", metavar="character"),
  make_option(c("-g", "--gtf"), type="character", help="GTF file", metavar="character")
)

opt <- parse_args(OptionParser(option_list=option_list))

# Function to run shell command and capture output
run_command <- function(cmd) {
  message(cmd)
  system(cmd)
}

path <- getwd()

# Read SRP file list
srp_list <- readLines(opt$srp_file)

##remove all MXE files in 

for (srp in srp_list) {
  message(srp)
  
  # Get rid of MXE, we aren't using them
  unlink(file.path(path, "rmats", opt$comparison, srp, "output", "*MXE*"))
  
  # Impute missing values
  message("impute missing")
  impute_cmd <- sprintf("Rscript scripts/impute.R %s %s %s output", srp, opt$comparison, path)
  run_command(impute_cmd)
  
  # Check if imputed file exists and run rmats_stat if it does
  imputed_file <- file.path(path, "rmats", opt$comparison, srp, "output", "imputed_JCEC.raw.input.SE.txt")
  if (file.exists(imputed_file)) {
    message("run rmats stat")
    
    for (event_type in c("A5SS", "A3SS", "SE", "RI")) {
      rmats_cmd <- sprintf(
        "python2.7 scripts/rMATS_unpaired.py %s %s 20 0.05 > %s_log.txt",
        file.path(path, "rmats", opt$comparison, srp, "output", sprintf("imputed_JCEC.raw.input.%s.txt", event_type)),
        file.path(path, "rmats", opt$comparison, srp, "output"),
        srp
      )
      run_command(rmats_cmd)
      file.rename(
        file.path(path, "rmats", opt$comparison, srp, "output", "rMATS_Result_P.txt"),
        file.path(path, "rmats", opt$comparison, srp, "output", sprintf("%s_rMATS_Result_P.txt", event_type))
      )
    }
    
    # Format files and fix p-values
    message("file formatting and fixing pvalues")
    format_cmd <- sprintf("Rscript scripts/reformat.R %s %s %s output", srp, opt$comparison, path)
    run_command(format_cmd)
  } else {
    message("imputed files do not exist")
  }
}
