
## Marta R. Hidalgo

library(hipathia)
context("Head function for matrices")

m <- matrix(0, ncol = nc, nrow = nr)

test_that("Performs well on dim(3,3)", {
    nr <- 3
    nc <- 3
    m <- matrix(0, ncol = nc, nrow = nr)
    hh <- hhead(m)
    expect_is(hh, "matrix")
    expect_equal(nrow(hh), nr)
    expect_equal(ncol(hh), nc)
    expect_true(all(hh == m[1:nrow(hh), 1:ncol(hh)]))
})

test_that("Performs well on dim(7,3)", {
    nr <- 7
    nc <- 3
    m <- matrix(0, ncol = nc, nrow = nr)
    hh <- hhead(m)
    expect_is(hh, "matrix")
    expect_equal(nrow(hh), nr)
    expect_equal(ncol(hh), nc)
    expect_true(all(hh == m[1:nrow(hh), 1:ncol(hh)]))
})

test_that("Performs well on dim(10,3)", {
    nr <- 10
    nc <- 3
    m <- matrix(0, ncol = nc, nrow = nr)
    hh <- hhead(m)
    expect_is(hh, "matrix")
    expect_equal(nrow(hh), 5)
    expect_equal(ncol(hh), nc)
    expect_true(all(hh == m[1:nrow(hh), 1:ncol(hh)]))
})

test_that("Performs well on dim(3,7)", {
    nr <- 3
    nc <- 7
    m <- matrix(0, ncol = nc, nrow = nr)
    hh <- hhead(m)
    expect_is(hh, "matrix")
    expect_equal(nrow(hh), nr)
    expect_equal(ncol(hh), nc)
    expect_true(all(hh == m[1:nrow(hh), 1:ncol(hh)]))
})

test_that("Performs well on dim(3,10)", {
    nr <- 3
    nc <- 10
    m <- matrix(0, ncol = nc, nrow = nr)
    hh <- hhead(m)
    expect_is(hh, "matrix")
    expect_equal(nrow(hh), nr)
    expect_equal(ncol(hh), 5)
    expect_true(all(hh == m[1:nrow(hh), 1:ncol(hh)]))
})

test_that("Performs well on dim(10,10)", {
    nr <- 10
    nc <- 10
    m <- matrix(0, ncol = nc, nrow = nr)
    hh <- hhead(m)
    expect_is(hh, "matrix")
    expect_equal(nrow(hh), 5)
    expect_equal(ncol(hh), 5)
    expect_true(all(hh == m[1:nrow(hh), 1:ncol(hh)]))
})
