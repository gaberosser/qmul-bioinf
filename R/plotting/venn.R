library(VennDiagram)

get_substr_by_index <- function(s, idx) {
  paste0(sapply(idx, FUN = function(n) substr(s, n, n)), collapse = '')
}


get_venn_area <- function(counts, n) {
  idx <- sapply(names(counts), FUN=function(x) get_substr_by_index(x, n))
  idx[idx == ""] <- "0"
  # all results in idx that do not contain any 0 characters
  do.call(sum, counts[grep("0", idx, invert = T)])
}


venn_diagram.from_blocks <- function(blocks, count_func=nrow, ...) {
  n = log2(length(blocks))
  if (n %% 1 != 0) {
    # TODO
    stop("Unexpected input length. We should have one entry for each possible Venn combination plus one called `contrasts`.")
  }
  if (n > 5) {
    stop(sprintf("Can't plot a Venn diagram with n>5 (supplied: %i)", n))
  }
  
  counts <- lapply(blocks, count_func)
  counts$contrasts <- NULL
  category = blocks$contrasts
  
  if (n >= 2) {
    area1 = get_venn_area(counts, 1)
    area2 = get_venn_area(counts, 2)
    n12 = get_venn_area(counts, c(1, 2))
    
    if (n == 2) {draw.pairwise.venn(area1 = area1, area2 = area2, cross.area = n12, category = category, ...)}
  }
  
  if (n >= 3) {
    area3 = get_venn_area(counts, 3)
    n13 = get_venn_area(counts, c(1,3))
    n23 = get_venn_area(counts, c(2,3))
    n123 = get_venn_area(counts, c(1,2,3))
    if (n == 3) {
      draw.triple.venn(
        area1 = area1, area2 = area2, area3 = area3,
        n12 = n12, n13 = n13, n23 = n23, n123 = n123,
        category = category, ...
      )
    }
    
  }
  
  if (n >= 4) {
    area4 = get_venn_area(counts, 4)
    n14 = get_venn_area(counts, c(1,4))
    n24 = get_venn_area(counts, c(2,4))
    n34 = get_venn_area(counts, c(3,4))
    n123 = get_venn_area(counts, c(1,2,3))
    n124 = get_venn_area(counts, c(1,2,4))
    n134 = get_venn_area(counts, c(1,3,4))
    n234 = get_venn_area(counts, c(2,3,4))
    n1234 = get_venn_area(counts, c(1,2,3,4))
    if (n == 4) {
      draw.quad.venn(area1=area1, area2=area2, area3=area3, area4=area4, 
                     n12=n12, n13=n13, n14=n14, n23=n23, n24=n24, n34=n34, 
                     n123=n123, n124=n124, n134=n134, n234=n234, 
                     n1234=n1234,
                     category=category, ...)
    }

  }
  
  if (n >= 5) {
    area5 = get_venn_area(counts, 5)
    n15 = get_venn_area(counts, c(1, 5))
    n25 = get_venn_area(counts, c(2, 5))
    n35 = get_venn_area(counts, c(3, 5))
    n45 = get_venn_area(counts, c(4, 5))
    n125 = get_venn_area(counts, c(1, 2, 5))
    n135 = get_venn_area(counts, c(1, 3, 5))
    n145 = get_venn_area(counts, c(1, 4, 5))
    n235 = get_venn_area(counts, c(2, 3, 5))
    n245 = get_venn_area(counts, c(2, 4, 5))
    n345 = get_venn_area(counts, c(3, 4, 5))
    n1235 = get_venn_area(counts, c(1, 2, 3, 5))
    n1245 = get_venn_area(counts, c(1, 2, 4, 5))
    n1345 = get_venn_area(counts, c(1, 3, 4, 5))
    n2345 = get_venn_area(counts, c(2, 3, 4, 5))
    n12345 = get_venn_area(counts, c(1, 2, 3, 4, 5))
    if (n == 5) {
      draw.quintuple.venn(
        area1, area2, area3, area4, area5, n12, n13, n14, n15,
        n23, n24, n25, n34, n35, n45, n123, n124, n125, n134,
        n135, n145, n234, n235, n245, n345, n1234, n1235,
        n1245, n1345, n2345, n12345,
        category = category, ...
      )
    }
  }
}
