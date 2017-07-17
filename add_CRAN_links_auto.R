path <- paste0(getwd(), "/README.Rmd")

rmd <- readLines(paste0(getwd(), "/README.Rmd"))

rmd2 <- gsub("\\* \\[(.+)]\\(\\)", "\\* \\[\\1]\\(https://cran.r-project.org/package=\\1\\)", rmd)

cat(paste(rmd2, collapse = "\n"), file = path)
