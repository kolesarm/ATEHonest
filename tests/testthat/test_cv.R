context("Check critical value calculations")

test_that("Limit casess for critical value", {
    expect_equal(cv(c(NA, 1, 100), alpha=0.1),
                 c(NA, cv(1, alpha=0.1), qnorm(1-0.1)+100))
    expect_equal(cv(0, alpha=0.2), qnorm(0.9))
})
