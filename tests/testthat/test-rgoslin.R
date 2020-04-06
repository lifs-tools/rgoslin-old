test_that("lipid name validation works", {
  expect_equal(rgoslin::isValidLipidName("PC(32:0)"), TRUE)
})

test_that("lipid name parsing works", {
  originalName <- "LPC(34:2;1)"
  vec <- rgoslin::parseLipidName(originalName)
  expect_equal(is.vector(vec), TRUE)
  expect_equal(vec[["Original Name"]], originalName)
  expect_equal(vec[["Normalized Name"]], "LPC 34:2;1")
  expect_equal(vec[["Lipid Maps Category"]], "GP")
  # expect_equal(vec[["Lipid Maps Category"]], "LPC 34:2;1")
})

test_that("multiple lipid names parsing works", {
  originalNames <- c("PC(32:0)", "LPC(34:2;1)", "TG(18:1_18:0_16:1)", "TAG 16:1/18:0/20:2")
  df <- rgoslin::parseLipidNames(originalNames)
  expect_equal(is.data.frame(df), TRUE)
  expect_equal(nrow(df), 4)
  expect_equal(as.character(df[1, "Original Name"]), originalNames[[1]])
  expect_equal(as.character(df[2, "Original Name"]), originalNames[[2]])
  expect_equal(as.character(df[3, "Original Name"]), originalNames[[3]])
  expect_equal(as.character(df[4, "Original Name"]), originalNames[[4]])
  expect_equal(as.character(df[1, "Normalized Name"]), "PC 32:0")
  expect_equal(as.character(df[2, "Normalized Name"]), "LPC 34:2;1")
  expect_equal(as.character(df[3, "Normalized Name"]), "TAG 18:1-18:0-16:1")
  expect_equal(as.character(df[4, "Normalized Name"]), "TAG 16:1/18:0/20:2")
})
