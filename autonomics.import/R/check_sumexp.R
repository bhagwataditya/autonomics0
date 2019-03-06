check_sumexp <- function(object){

   i <- 0
   msg <- ''

   # subgroup var
   #-------------
   if (!'subgroup' %in% svars(object)){
      i <- i+1
      msg <- sprintf("%s\n   %d. Add 'subgroup'    : object$subgroup <- subgroup_values", msg, i)

      # subgroup values
      #----------------
   } else if (any(is.na(sdata(object)$subgroup))){
      i <- i+1
      msg <- sprintf("%s\n   %d. Complete 'subgroup': object$subgroup contains missing values", msg, i)
   }

   # block var
   #----------
   if (!'block' %in% svars(object)){
      i <- i+1
      msg <- sprintf("%s\n   %d. Add 'block'       : object$block <- block_values (if data has repeated measures)", msg, i)

      # block values
      #-------------
   } else if (any(is.na(sdata(object)$block))){
      i <- i+1
      msg <- sprintf("%s\n   %d. Complete 'block': object$block contains missing values", msg, i)
   }

   if (i>0){
      autonomics.support::cmessage(sprintf('Action required%s', msg))
   }

}
