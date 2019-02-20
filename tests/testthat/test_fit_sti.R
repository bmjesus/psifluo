context("testing the sti_fit output")
library(psifluo)
test_that("fit_sti output is ok", {
  out_sti <- fit_sti(x = 0.0002167971, y = sti_data$light_1,fit_model = 'all')
  expect_equal(out_sti[[1]], 0.0002167971)
  expect_is(out_sti[[2]], "list")

})



