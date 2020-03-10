test_that("lipid name validation works", {
  expect_equal(rgoslin::isValidLipidName("PC(32:0)"), TRUE)
})

test_that("lipid name parsing works", {
  originalName <- "LPC(18:0;1_16:2)"
  vec <- rgoslin::parseLipidName(originalName)
  expect_equal(is.vector(vec), TRUE)
  expect_equal(vec[["Original Name"]], originalName)
  expect_equal(vec[["Normalized Name"]], "LPC 18:0;1_16:2")
  expect_equal(vec[["Lipid Maps Category"]], "GP")
  # expect_equal(vec[["Lipid Maps Category"]], "LPC 34:2;1")
})
