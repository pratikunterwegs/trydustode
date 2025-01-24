test_that("dust system works", {
  sys <- dust2::dust_system_create(sirode, list(1, 1e6, 0.3, 0.2, 0.1))
  dust2::dust_system_set_state_initial(sys)

  times <- 1:100
  state <- dust2::dust_system_simulate(sys, times)
  state <- dust2::dust_unpack_state(sys, state)

  expect_lte(
    max(diff(state$S)), 0.0
  )
})

test_that("using dust in function", {
  expect_no_condition(
    f()
  )
})
