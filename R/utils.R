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
  # Melt the data.table so all values are in one column called "value".
  dat <- melt(dat, id.vars = "agg_var")
  # Cast the data.table back into the original shape, and aggregate.
  dat <- dcast.data.table(
    dat, agg_var ~ variable, value.var = "value",
    fun.aggregate = median, na.rm = TRUE
  )
  rownames(dat) <- dat$agg_var
  return(dat)
}

mean_by <- function(dat, xs) {
  # Convert to data.table.
  dat <- data.table(dat)
  # Append the vector of group names as an extra column.
  dat$agg_var <- xs
  # Melt the data.table so all values are in one column called "value".
  dat <- melt(dat, id.vars = "agg_var")
  # Cast the data.table back into the original shape, and aggregate.
  dat <- dcast.data.table(
    dat, agg_var ~ variable, value.var = "value",
    fun.aggregate = mean, na.rm = TRUE
  )
  rownames(dat) <- dat$agg_var
  return(dat)
}

min_by <- function(dat, xs) {
  # Convert to data.table.
  dat <- data.table(dat)
  # Append the vector of group names as an extra column.
  dat$agg_var <- xs
  # Melt the data.table so all values are in one column called "value".
  dat <- melt(dat, id.vars = "agg_var")
  # Cast the data.table back into the original shape, and aggregate.
  dat <- dcast.data.table(
    dat, agg_var ~ variable, value.var = "value",
    fun.aggregate = min, na.rm = TRUE
  )
  rownames(dat) <- dat$agg_var
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
  dat <- subset(dat, select = -c(agg_var))
  return(dat)
}