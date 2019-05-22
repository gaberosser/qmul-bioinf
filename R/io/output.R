source('_settings.R')

# create output directory if it doesn't exist
dir.create(out.dir, showWarnings = FALSE)


getOutputDir <- function(name) {
  today <- Sys.Date()
  date_str <- format(today, format="%Y-%m-%d")
  
  name.num <- paste(name, date_str, '{n}', sep = '.')

  i <- 0
  ff <- file.path(out.dir, sub('{n}', i, name.num, fixed=T))
  while (file.exists(ff)) {
    i <- i + 1
    ff <- file.path(out.dir, sub('{n}', i, name.num, fixed=T))
  }
  dir.create(ff)
  message(paste0("Created output directory ", ff))
  return(ff)
}
