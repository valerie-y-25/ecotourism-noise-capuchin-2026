# ICC TEST: Here is where we tested our Practice BORIS Data. Raw Practice BORIS Data is in the "BORIS Focal Data" Folder in BOX 
ratings <- data.frame(
  rater1 = c(10.167, 24.637, 8.298, 7.385),
  rater2 = c(10.382, 24.897, 10.146, 7.706),
  rater3 = c(8.311, 26.252, 11.994, 6.749)
)
icc_result <- ICC(ratings)
print(icc_result)