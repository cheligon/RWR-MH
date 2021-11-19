get.data.from.zip.at.url <- function(zipfile.url, dir, file.sep = "\t"){
  zip.path <- paste(dir, basename(zipfile.url), sep = "/")
  if (!file.exists(zip.path)){
    download.file(zipfile.url, zip.path, method = "curl")  
  }
  internal.file <- tools::file_path_sans_ext(basename(zip.path))
  unzipped.file <- unz(zip.path, internal.file)
  expression.table <- read.table(unzipped.file, header = TRUE, sep = file.sep)
  return(expression.table)
}