test_that("prepare_actt_data works", {
  data("dummy_actt_ch")
  data("dummy_actt_long")

  result <- prepare_actt_data()
  expect_true(is.list(result))
  expect_true("data.w" %in% names(result))
  expect_s3_class(result$data.w, "data.frame")
  expect_true(nrow(result$data.w) > 0)
})
