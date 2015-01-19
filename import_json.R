import_json <- function (json_file) {
  library('rjson')
  json_file <- fromJSON(file = json_file, method='C')
  unlisted <- unname(unlist(json_file))
  
  df<- (as.data.frame(unlisted))
  sequence1<- seq(1,nrow(df),2)
  sequence2<- seq(2,nrow(df),2)
  reps <- as.data.frame (df[sequence1,])
  conds<- as.data.frame (df[sequence2,])
  output_frame <- cbind(reps,conds) 
  print(output_frame)
}