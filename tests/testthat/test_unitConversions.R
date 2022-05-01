#require(testthat)
context("unitConversions")


test_that("e.to.rH",{
  Tair <- 25
  eSat <- Esat.slope(Tair)$Esat
  e <- seq(0,eSat, length.out = 5)
  rH0 <- VPD.to.rH(e.to.VPD(e, Tair),Tair)
  rH <- e.to.rH(e, Tair)
  expect_equal(rH, rH0)
  #
  expect_warning(
    rHOversat <- e.to.rH(eSat + 1e-3, Tair)
  )
  expect_equal(rHOversat, 1)
})

test_that("kg.to.mol",{
  mass <- 10
  molarMass <- bigleaf.constants()$H2Omol
  amountOfSubstance <- kg.to.mol(mass, molarMass)
  expect_equal(amountOfSubstance, mass/molarMass)
})

