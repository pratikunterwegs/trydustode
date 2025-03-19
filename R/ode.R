#' Small SIR function using dust2
#'
#' @param times A vector of times.
#' @export
fsir <- function(times = 0:10) {
  sys <- dust2::dust_system_create(sirode)
  dust2::dust_system_set_state_initial(sys)
  state <- dust2::dust_system_simulate(sys, times)
  state <- dust2::dust_unpack_state(sys, state)

  state
}

#' Bouncing ball simulation
#'
#' @param times A vector of times.
#' @export
fball <- function(times = 0:10, h = 100) {
  sys <- dust2::dust_system_create(bounce, pars = list(height = h))

  dust2::dust_system_set_state_initial(sys)
  state <- dust2::dust_system_simulate(sys, times)
  state <- dust2::dust_unpack_state(sys, state)

  state
}
