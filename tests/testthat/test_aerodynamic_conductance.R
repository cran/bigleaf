context("aerodynamic conductance")


test_that("aerodynamic conductance",{
  df       <- data.frame(Tair=25,pressure=100,wind=c(3,4,5),ustar=c(0.5,0.6,0.65),H=c(200,230,250)) 
  
  Ga_Thom  <- aerodynamic.conductance(df,Rb_model="Thom_1972")[,c("Ga_m","Ga_h","kB_h")]
  Ga_Thom2 <- aerodynamic.conductance(df,Rb_model="Thom_1972",zr=40,zh=25,d=17.5,z0m=2,wind_profile=TRUE)[,c("Ga_m","Ga_h","kB_h")]
  Ga_Su    <- aerodynamic.conductance(df,Rb_model="Su_2001",zr=40,zh=25,d=17.5,Dl=0.05,N=2,fc=0.8)[,c("Ga_m","Ga_h","kB_h")]
  
  exp_Ga_Thom  <- data.frame("Ga_m"=c(0.0833333333333,0.0900000000000,0.0845000000000),
                             "Ga_h"=c(0.0457788244529,0.0504335618459,0.0497559607901),
                             "kB_h"=c(2.01805295243,2.14437093732,2.20229602024))
  exp_Ga_Thom2 <- data.frame("Ga_m"=c(0.177067652414,0.175231551534,0.178738919785),
                             "Ga_h"=c(0.0645506129505,0.0693303869614,0.0721577391297),
                             "kB_h"=c(2.01805295243,2.14437093732,2.20229602024))
  exp_Ga_Su    <- data.frame("Ga_m"=c(0.0833333333333,0.0900000000000,0.0845000000000),
                             "Ga_h"=c(0.0412084575303,0.0426710767005,0.0415922922898),
                             "kB_h"=c(2.51470694818,3.03169579229,3.25359102569))
  
  expect_that(Ga_Thom,equals(exp_Ga_Thom,tolerance=1e-9))
  expect_that(Ga_Thom2,equals(exp_Ga_Thom2,tolerance=1e-9))
  expect_that(Ga_Su,equals(exp_Ga_Su,tolerance=1e-9))
})


test_that("stability correction",{
  zeta <- c(-1.5,-1.0,NA,0,0.2,0.5)
  
  zeta1 <- stability.correction(zeta,formulation="Dyer_1970")
  zeta2 <- stability.correction(zeta,formulation="Businger_1971")
  
  exp_zeta1 <- data.frame("psi_h"=c(2.19722457734,1.88122728421,NA,0.00000000000,-1.00000000000,-2.50000000000),
                          "psi_m"=c(2.06103593879,1.77180302576,NA,0.00000000000,-1.00000000000,-2.50000000000))
  exp_zeta2 <- data.frame("psi_h"=c(1.86237682147,1.56422247632,NA,0.00000000000,-1.56000000000,-3.90000000000),
                          "psi_m"=c(2.19971077502,1.90366580737,NA,0.00000000000,-1.20000000000,-3.00000000000))
  
  expect_that(zeta1,equals(exp_zeta1,tolerance=1e-9))
  expect_that(zeta2,equals(exp_zeta2,tolerance=1e-9))
    
})