library(dplyr)

# Check if arguments are passed
if (length(commandArgs(trailingOnly = TRUE)) < 3) {
  stop("Insufficient arguments. Include SRP# and library type (single or paired) and comparison: Rscript script.R SRPXXXXXX single AD_v_healthy")
}

args <- commandArgs(trailingOnly = TRUE)

srp <- args[1]
library <- args[2]
comp <- args[3]

dir.create("rmats")
dir.create("rmats/tmp")
dir.create(paste0("rmats/",comp))
dir.create(paste0("rmats/",comp, "/", srp))
dir.create(paste0("rmats/",comp, "/", srp, "/tmp"))
dir.create(paste0("rmats/",comp, "/", srp, "/input"))
dir.create(paste0("rmats/",comp, "/", srp, "/output"))


##write input and output .txt files containing bam paths
metadata <- readRDS("metadata.rds")
metadata <- metadata %>% filter(SRA.STUDY == srp)
metadata$path <- paste0("/data/", srp, "/bam/", metadata$RUN, ".bam")

metadata_control <- metadata %>% filter(GROUP == "group2")
metadata_condition <- metadata %>% filter(GROUP == "group1")

control_bam <- metadata_control$path
control_bam <- paste(control_bam, collapse = ",")

condition_bam <- metadata_condition$path
condition_bam <- paste(condition_bam, collapse = ",")

write(control_bam, file = paste0("rmats/",comp, "/", srp, "/input/control.txt"))
write(condition_bam, file = paste0("rmats/",comp, "/", srp, "/input/condition.txt"))

#get read length
read_length_path <- paste0("/data/", srp, "/multiqc/multiqc_data/multiqc_fastqc.txt")
multiqc_fastqc <- read.delim(read_length_path)

read_length <- round(sum(multiqc_fastqc$median_sequence_length)/nrow(multiqc_fastqc))


rmats_command_prep <- paste(
  "rmats.py",
  "--b1", paste0("rmats/",comp, "/", srp, "/input/condition.txt"),
  "--b2", paste0("rmats/",comp, "/", srp, "/input/control.txt"),
  "--gtf /data/Homo_sapiens.GRCh38.103.gtf ",
  "--od", paste0("rmats/",comp, "/", srp, "/output/"),
  "-t", library,
  "--readLength", read_length,
  "--variable-read-length",
  "--allow-clipping",
  "--tmp", paste0("rmats/",comp, "/", srp, "/tmp/"),
  "--nthread 20",
  "--task prep"
)

rmats_command_post <- paste(
  "rmats.py",
  "--b1", paste0("rmats/",comp, "/", srp, "/input/condition.txt"),
  "--b2", paste0("rmats/",comp, "/", srp, "/input/control.txt"),
  "--gtf /data/Homo_sapiens.GRCh38.103.gtf ",
  "--od", paste0("rmats/",comp, "/", srp, "/output/"),
  "-t", library,
  "--readLength", read_length,
  "--variable-read-length",
  "--allow-clipping",
  "--tmp", paste0("rmats/",comp, "/", srp, "/tmp/"),
  "--nthread 20",
  "--task post"
)

print(rmats_command_prep)
system(rmats_command_prep)


print(rmats_command_post)
system(rmats_command_post)