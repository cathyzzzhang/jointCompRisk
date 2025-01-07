test_that("prep_data_wrt_cif works", {
  res_part1 <- prep_data_cif()
  res_part2 <- prep_data_wrt_cif(part1_output = res_part1)
  expect_true(is.list(res_part2))
  expect_true("data.ws.death" %in% names(res_part2))
  expect_s3_class(res_part2$data.ws.death, "data.frame")
  expect_true(nrow(res_part2$data.ws.death) > 0)
})
