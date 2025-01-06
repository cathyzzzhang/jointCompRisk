test_that("w_surv basic usage", {
  data <- data.frame(
    D_time = c(5, 10, 15),
    D_status = c(1, 0, 1),
    etype = c(1, 2, 1),
    wU = c(10, 15, 20)
  )

  result <- w_surv(data, t = 15, eta = 1, m = 20)
  expect_true(is.numeric(result))
  expect_gte(result, 0)
})
