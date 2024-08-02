assertNames <- function(obj1, obj2) {
  if (!all(names(obj1) %in% names(obj2)) ||
      !all(names(obj2) %in% names(obj1)) {
        stop("Object `names()` do not match")
      }
}
