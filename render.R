# From Tyler Rinker: https://github.com/trinker/testing_Rmd

rmarkdown::render("README.Rmd", "all")

md_toc <- function(path = paste0(getwd(), "/README.md")){
  x <- suppressWarnings(readLines(path))
  inds <- 1:(which(!grepl("^\\s*-", x))[1] - 1)
  temp <- gsub("(^[ -]+)(.+)", "\\1", x[inds])
  content <- gsub("^[ -]+", "", x[inds])
  x[inds] <- sprintf("%s[%s](#%s)", temp, content, gsub("\\s", "-", tolower(content)))
  cat(paste(x, collapse = "\n"), file = path)
}

md_toc()
