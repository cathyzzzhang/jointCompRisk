test_that("prep_data_cif works", {
  res_part1 <- prep_data_cif()
  expect_true(is.list(res_part1))
  expect_true("data.w" %in% names(res_part1))
  expect_s3_class(res_part1$data.w, "data.frame")
  expect_true(nrow(res_part1$data.w) > 0)
})
