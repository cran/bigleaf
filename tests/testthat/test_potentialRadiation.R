#require(testthat)
context("potentialRadiation")


test_that("Regression test potential.radiation",{
  hour <- seq(8,16, by = 1)
  potRadSolar <- potential.radiation(160, hour, 39.94, -5.77, timezone = +1)
  potRadLocal <- potential.radiation(
    160, hour, 39.94, -5.77, timezone = +1, useSolartime = FALSE)
  expSolar <- structure(c(
    484.152670743821, 717.876981534078, 925.130678985721,
    1091.78976612035, 1206.4967015669, 1261.43439723686, 1252.85893995917,
    1181.35473297567, 1051.79466982602))
  expLocal <- structure(c(
    797.589402859243, 991.498827921097, 1140.29076299454,
    1233.82528359462, 1265.72816671554, 1233.82528359462, 1140.29076299454,
    991.498827921097, 797.589402859243))
  expect_that( potRadSolar, equals(expSolar) )
  expect_that( potRadLocal, equals(expLocal) )
})

test_that("potential.radiation: warn on non-matching day-hour length",{
  hour <- seq(8,16, by = 0.1)
  expect_warning(
    potRadSolar <- potential.radiation(160:161, hour, 39.94, -5.77, timezone = +1)
  )
})

