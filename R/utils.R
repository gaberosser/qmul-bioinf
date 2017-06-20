#' Bind columns, adding additional rows if there is any mismatch.
#' Missing values are filled with NA.
cbind.outer <- function(...) {
  dat <- NULL
  
  for (t in list(...)) {
    if (is.null(dat)) {
      dat <- data.frame(t)
    } else {
      this <- data.frame(t)
      if (
        (length(rownames(t)) != length(rownames(dat))) | (length(intersect(rownames(t), rownames(dat))) != length(t))
      ) {
        # expand the list of rows
        new_rows <- setdiff(rownames(this), rownames(dat))
        dat[new_rows,] <- NA
        new_rows <- setdiff(rownames(dat), rownames(this))
        this[new_rows,] <- NA
      }
      dat <- cbind(dat, this)
    }
  }
  dat
}

#' Bind rows, adding additional columns if there is any mismatch.
#' Missing values are filled with NA.
rbind.outer <- function(...) {
  dat <- NULL
  
  for (t in list(...)) {
    if (is.null(dat)) {
      dat <- data.frame(t)
    } else {
      this <- data.frame(t)
      if (
        (length(colnames(t)) != length(colnames(dat))) | (length(intersect(colnames(t), colnames(dat))) != length(t))
      ) {
        # expand the list of rows
        new_rows <- setdiff(colnames(this), colnames(dat))
        dat[,new_rows] <- NA
        new_rows <- setdiff(colnames(dat), colnames(this))
        this[,new_rows] <- NA
      }
      dat <-rbind(dat, this)
    }
  }
  dat
}

#' Replicate a vector n times and stack as columns to form a matrix
rep.col<-function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

#' Replicate a vector n times and stack as rows to form a matrix
rep.row<-function(x,n){
  matrix(rep(x,each=n),nrow=n)
}

#' Function source: http://slowkow.com/notes/data-table-aggregate/
#' Author: Kamil Slowikowski
#' Aggregate all columns of a matrix or dataframe, where rows are
#' aggregated by a vector of values. 100 times faster than stats::aggregate.
#'
#' @param dat A numeric matrix or data.frame.
#' @param xs A vector of groups (e.g. gene names).
#' @return A data.table with the aggregated mean for each group.
#' @seealso stats::aggregate
median_by <- function(dat, xs) {
  # Convert to data.table.
  dat <- data.table(dat)
  # Append the vector of group names as an extra column.
  dat$agg_var <- xs
  # Remove null entries in the aggregation variable, as these cannnot be used for indexing
  dat <- dat[!is.na(dat$agg_var),]  
  # Melt the data.table so all values are in one column called "value".
  dat <- melt(dat, id.vars = "agg_var")
  # Cast the data.table back into the original shape, and aggregate.
  dat <- dcast.data.table(
    dat, agg_var ~ variable, value.var = "value",
    fun.aggregate = median, na.rm = TRUE
  )
  # set rownames by agg_var and remove the unneeded column
  rownames(dat) <- dat$agg_var
  dat[,agg_var:=NULL]
  return(dat)
}

mean_by <- function(dat, xs) {
  # Convert to data.table.
  dat <- data.table(dat)
  # Append the vector of group names as an extra column.
  dat$agg_var <- xs
  # Remove null entries in the aggregation variable, as these cannnot be used for indexing
  dat <- dat[!is.na(dat$agg_var),]  
  # Melt the data.table so all values are in one column called "value".
  dat <- melt(dat, id.vars = "agg_var")
  # Cast the data.table back into the original shape, and aggregate.
  dat <- dcast.data.table(
    dat, agg_var ~ variable, value.var = "value",
    fun.aggregate = mean, na.rm = TRUE
  )
  # set rownames by agg_var and remove the unneeded column
  rownames(dat) <- dat$agg_var
  dat[,agg_var:=NULL]
  return(dat)
}

min_by <- function(dat, xs) {
  # Convert to data.table.
  dat <- data.table(dat)
  # Append the vector of group names as an extra column.
  dat$agg_var <- xs
  # Remove null entries in the aggregation variable, as these cannnot be used for indexing
  dat <- dat[!is.na(dat$agg_var),]  
  # Melt the data.table so all values are in one column called "value".
  dat <- melt(dat, id.vars = "agg_var")
  # Cast the data.table back into the original shape, and aggregate.
  dat <- dcast.data.table(
    dat, agg_var ~ variable, value.var = "value",
    fun.aggregate = min, na.rm = TRUE
  )
  # set rownames by agg_var and remove the unneeded column
  rownames(dat) <- dat$agg_var
  dat[,agg_var:=NULL]
  return(dat)
}

max_by <- function(dat, xs) {
  # Convert to data.table.
  dat <- data.table(dat)
  # Append the vector of group names as an extra column.
  dat$agg_var <- xs
  # Remove null entries in the aggregation variable, as these cannnot be used for indexing
  dat <- dat[!is.na(dat$agg_var),]
  # Melt the data.table so all values are in one column called "value".
  dat <- melt(dat, id.vars = "agg_var")
  # Cast the data.table back into the original shape, and aggregate.
  dat <- dcast.data.table(
    dat, agg_var ~ variable, value.var = "value",
    fun.aggregate = max, na.rm = TRUE
  )
  # set rownames by agg_var and remove the unneeded column
  rownames(dat) <- dat$agg_var
  dat[,agg_var:=NULL]
  return(dat)
}

#' Given a variable number of vector inputs, generate a list with names giving the boolean intersection of the inputs
#' For example, out['101'] contains entries that are present in the first AND third input, but NOT the second
#' The inputs are vectors or factors (?) that can be passed through intersection(). 
venn_sets <- function(...) {
  ids <- list(...)
  blocks <- list()
  
  nbit <- length(ids)
  n <- 2 ** nbit - 1
  
  
  for (i in seq(1, n)) {
    comb <- as.integer(intToBits(i))[1:nbit]
    idx.in <- which(comb == 1)
    idx.out <- which(comb == 0)
    
    ids.in <- Reduce(intersect, lapply(idx.in, function(x){ ids[[x]] }))
    ids.out <- Reduce(union, lapply(idx.out, function(x){ ids[[x]] }))
    ids.this <- setdiff(ids.in, ids.out)

    blocks[[paste0(comb, collapse = '')]] <- ids.this
  }
  
  blocks
}