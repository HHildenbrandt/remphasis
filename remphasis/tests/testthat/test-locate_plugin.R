context("locate_plugin")

testthat::test_that("test_locate_plugin", {
  testthat::skip_on_os("solaris") # solaris is often very bitchy, better avoid it.
  
  testthat::expect_warning(locate_plugin("noname"))

 # uncomment code once the library below is available (e.g. on a non-private github)
 # library(remphasis_rpd1)
 # testthat::expect_silent(locate_plugin("rpd1"))
    
})
