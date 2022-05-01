context("meteorological variables")

test_that("wet bulb temperature",{
  Tair  <- seq(5,45,5)
  Tw    <- wetbulb.temp(Tair=Tair,pressure=100,VPD=0.5,accuracy=1e-06,Esat.formula="Sonntag_1990")
  expTw <- c(0.590813, 6.273059, 11.898636, 17.448358, 22.915922, 
             28.304377, 33.622011, 38.879140, 44.086124)
  expect_that(Tw,equals(expTw) )
})