test_that("lipid name validation works", {
  expect_equal(rgoslin::isValidLipidName("PC(32:0)"), TRUE)
})

test_that("lipid name parsing works", {
  originalName <- "LPC(18:0;1_16:2)"
  vec <- rgoslin::parseLipidName(originalName)
  expect_equal(is.vector(vec), TRUE)
  expect_equal(vec[["originalName"]], originalName)
  expect_equal(vec[["nativeLevelName"]], "LPC 18:0;1_16:2")
  expect_equal(vec[["category"]], "GP")
  expect_equal(vec[["species"]], "LPC 34:2;1")
})
