test_that("sample rules drop always-treated units and set the period step", {
  id <- rep(1:4, each = 4)
  t <- rep(c(2, 4, 6, 8), 4)
  g <- rep(c(0, 2, 6, 6), each = 4)
  timing <- jw_timing(id, t, g, anticipation = 1L)
  expect_equal(timing$keep, g != 2)
  expect_equal(timing$gap, 2)
  expect_equal(timing$gmax, Inf)
})

test_that("without never-treated units the last cohort becomes the control", {
  id <- rep(1:3, each = 4)
  t <- rep(1:4, 3)
  g <- rep(c(2, 3, 4), each = 4)
  timing <- jw_timing(id, t, g, anticipation = 1L)
  expect_equal(timing$gmax, 4)
  expect_equal(timing$keep, t < 4)
})

test_that("anticipation removes units first seen inside the anticipation window", {
  id <- rep(1:3, c(4, 4, 3))
  t <- c(1:4, 1:4, 2:4)
  g <- rep(c(0, 3, 3), c(4, 4, 3))
  expect_equal(jw_timing(id, t, g, anticipation = 2L)$keep, id != 3)
})

test_that("cells follow the treatment definitions in cohort-period order", {
  t <- rep(1:4, 3)
  g <- rep(c(0, 3, 4), each = 4)
  never <- jw_cells(t, g, gap = 1, gmax = Inf, anticipation = 1L, never = TRUE)
  expect_equal(levels(never$.jw_cell), c("base", "g3_t1", "g3_t3", "g3_t4", "g4_t1", "g4_t2", "g4_t4"))
  notyet <- jw_cells(t, g, gap = 1, gmax = Inf, anticipation = 1L, never = FALSE)
  expect_equal(levels(notyet$.jw_cell), c("base", "g3_t3", "g3_t4", "g4_t4"))
  expect_equal(notyet$.jw_etr, g > 0 & t >= g)
})

test_that("cells sort periods numerically", {
  t <- rep(c(1, 4, 7, 10), 2)
  g <- rep(c(0, 7), each = 4)
  cells <- jw_cells(t, g, gap = 3, gmax = Inf, anticipation = 1L, never = TRUE)
  expect_equal(levels(cells$.jw_cell), c("base", "g7_t1", "g7_t7", "g7_t10"))
})

test_that("the anticipation window follows the period step on irregular grids", {
  t <- rep(c(2, 5, 7, 9), 2)
  g <- rep(c(0, 5), each = 4)
  notyet <- jw_cells(t, g, gap = 2, gmax = Inf, anticipation = 2L, never = FALSE)
  expect_equal(levels(notyet$.jw_cell), c("base", "g5_t5", "g5_t7", "g5_t9"))
})
