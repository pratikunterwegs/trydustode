## Generated by dust2 (version 0.3.16) - do not edit
sirode <- structure(
  function() get("sirode"),
  class = "dust_system_generator",
  name = "sirode",
  package = "trydustode",
  path = NULL,
  parameters = data.frame(
    name = c("I0", "N", "beta", "sigma", "gamma", "n_strata"),
    type = c("real_type", "real_type", "real_type", "real_type", "real_type", "int"),
    constant = c(FALSE, TRUE, FALSE, FALSE, FALSE, FALSE)),
  properties = list(
    time_type = "continuous",
    has_compare = TRUE,
    has_adjoint = FALSE),
  default_dt = NULL)
