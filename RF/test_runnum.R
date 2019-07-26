print("script running")

runnum  <- Sys.getenv(c('runnum'))
print(runnum)
a=runnum
write.table(a, paste0("project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/test/",runnum,".txt"))