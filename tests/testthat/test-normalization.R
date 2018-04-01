
## Marta R. Hidalgo

library(hipathia)
context("Normalize_data function")

data("brca_data")
exp_data <- normalize_data(brca_data)

test_that("Resulting object is a matrix", {
    expect_is(exp_data, "matrix")
})

test_that("Exp.data parameter is necessary", {
    expect_error(normalize_data())
})

test_that("Truncation.percentil must be in [0,1]", {
    expect_error(normalize_data(brca_data, truncation_percentil = 5))
    expect_error(normalize_data(brca_data, truncation_percentil = -5))
    expect_error(normalize_data(brca_data, truncation_percentil = -0.1))
})

test_that("Dimnames are preserved", {
    expect_equal(colnames(brca_data), colnames(exp_data))
    expect_equal(rownames(brca_data), rownames(exp_data))
})

test_that("All values in [0,1]", {
    expect_true(all(exp_data >= 0 & exp_data <= 1))
    expect_equal(0, sum(is.na(exp_data)))
})

test_that("All distributions are the same when by.quantiles is TRUE", {
    exp_data <- normalize_data(brca_data, by_quantiles = T)
    q <- quantile(exp_data[,1])
    expect_true(all(apply(exp_data, 2, function(x) all(quantile(x) == q))))
})
