# General code for catching an error

# simple example that will generate an error

myfunction <- function(...){

  for(i in 1:10){
    if(i < 5){
      print(i)
    }else{
      blabla # will cause error: object blabla not found
    }
  }
 return('ok')
}

# wrong way: without trycatch. This will not save the output on line 26.
#out <- myfunction()

# with try catch, this will actually save the output on line 26
out <- tryCatch(expr = myfunction(),
                error = function(e){conditionMessage(e)})



save(out, file = '~/Downloads/trycatch.Rdata')
