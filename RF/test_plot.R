
x = c(1, 4, 5, 10, 22, 24)
y = c(1, 200, 400, 500, 550, 600)

plot1 = plot(x, y)

jpeg("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_2/testplot.jpeg", width = 500, height = 500)
plot(x, y)
dev.off()

save.image(file = "/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_2/test_plot.RData")
