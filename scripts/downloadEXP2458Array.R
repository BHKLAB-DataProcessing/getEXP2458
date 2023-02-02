require(downloader)

args <- commandArgs(trailingOnly = TRUE)
output_dir <- args[1]

zipdir1 <- "https://www.ebi.ac.uk/arrayexpress/files/E-MEXP-2458/E-MEXP-2458.raw.1.zip"
zipdir2 <- "https://www.ebi.ac.uk/arrayexpress/files/E-MEXP-2458/E-MEXP-2458.raw.2.zip"
dir3 <- "https://www.ebi.ac.uk/arrayexpress/files/E-MEXP-2458/E-MEXP-2458.sdrf.txt"

# require(R.utils) || stop("Library R.utils is not available!")

dwl.status1 <- download(zipdir1, destfile = file.path(output_dir, "E-MEXP-2458.raw.1.zip"))
dwl.status2 <- download(zipdir2, destfile = file.path(output_dir, "E-MEXP-2458.raw.2.zip"))
dwl.status3 <- download(dir3, destfile = file.path(output_dir, "E-MEXP-2458.sdrf.txt"))

unzip(file.path(output_dir, "E-MEXP-2458.raw.1.zip"), exdir = file.path(output_dir, "raw"))
unzip(file.path(output_dir, "E-MEXP-2458.raw.2.zip"), exdir = file.path(output_dir, "raw"))

unlink(file.path(output_dir, "E-MEXP-2458.raw.1.zip"))
unlink(file.path(output_dir, "E-MEXP-2458.raw.2.zip"))
