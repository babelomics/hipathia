
## Marta R. Hidalgo

library(hipathia)
context("Wilcoxon")

data("path_vals")
data("brca")
sample_group <- colData(brca)[,1]

test_that("Resulting object is a data.frame", {
    expect_is(comp, "data.frame")
})

test_that("Output is correct", {
    comp <- do_wilcoxon(path_vals, sample_group, g1 = "Tumor", g2 = "Normal")
    expect_equal(5, ncol(comp))
    expect_equal(nrow(comp), nrow(path_vals))
    expect_equal(rownames(comp), rownames(path_vals))
    expect_equal("numeric", class(comp$p.value))
    expect_equal("numeric", class(comp$statistic))
    expect_equal("numeric", class(comp$FDRp.value))
    expect_true(all(comp$`UP/DOWN` %in% c("UP", "DOWN")))
})

test_that("Equal datasets are not significant", {
    rnormdata <- rnorm(n = nrow(path_vals) * 20)
    eq_vals1 <- matrix(rnormdata,
                       nrow = nrow(path_vals),
                       ncol = 20,
                       dimnames = list(rownames(path_vals),
                                       colnames(path_vals)[1:20]))
    eq_vals2 <- matrix(rnormdata,
                       nrow = nrow(path_vals),
                       ncol = 20,
                       dimnames = list(rownames(path_vals),
                                       colnames(path_vals)[21:40]))
    eq_vals <- cbind(eq_vals1, eq_vals2)
    comp <- do.wilcoxon(eq_vals, sample_group, g1 = "Tumor", g2 = "Normal")

    expect_true(all(comp$p.value == 1))
    expect_true(all(comp$FDRp.value == 1))
    expect_true(all(comp$statistic == 0))
})

