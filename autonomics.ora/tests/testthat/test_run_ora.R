library(testthat)

context("run_ora")

#-------------------------
#         Pathway1
#          Y   N
#        ---------------
# Query Y| 100 100 | 200
#       N|  0  800 | 800
#         --------------
#         100  900 | 1000
#------------------------
query <- sprintf('s%02d', 0:199)
universe <- c(query, sprintf('s%02d', 200:999))
pathway_list <- list(pathway1 = sprintf('s%02d', 0:99))
test_that('run_ora returns expected results for test example 1', {
  a <- run_ora(query, universe, pathway_list)$p
  b <- phyper(100-1, 100, 900, 200, lower.tail = FALSE)
  c <- fisher.test(matrix(c(100, 0, 100, 800), nrow = 2))$p.value
  expect_equal(a,b)
  expect_equal(a,c)
})


#-------------------------
#         Pathway2
#          Y   N
#        ---------------
# Query Y| 50  150 | 200
#       N| 50  750 | 800
#         --------------
#         100  900 | 1000
#------------------------
query <- sprintf('s%02d', 0:199)
universe <- c(query, sprintf('s%02d', 200:999))
pathway_list <- list(pathway2 = sprintf('s%02d', 150:249))
test_that('run_ora returns expected results for test example 2', {
  a <- run_ora(query, universe, pathway_list)$p
  b <- phyper(50-1, 100, 900, 200, lower.tail = FALSE)
  c <- fisher.test(matrix(c(50, 50, 150, 750), nrow = 2))$p.value
  expect_equal(a,b)
  expect_equal(a,c)
})

  