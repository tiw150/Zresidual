
print.Zresidual <- function(x){
  UseMethod("print")
  print(as.vector(x))
}
