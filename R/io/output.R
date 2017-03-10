source('_settings.R')

# create output directory if it doesn't exist
dir.create(out.dir, showWarnings = FALSE)


getOutputDir <- function(name) {
  name.num <- '{base}.{n}'
  name.num <- sub('{base}', name, name.num, fixed=T)
  i <- 0
  ff <- file.path(out.dir, sub('{n}', i, name.num, fixed=T))
  while (file.exists(ff)) {
    i <- i + 1
    ff <- file.path(out.dir, sub('{n}', i, name.num, fixed=T))
  }
  dir.create(ff)
  print(paste0("Created output directory ", ff))
  return(ff)
}
