#include <cpp11.hpp>
#include <cpp11eigen.hpp>
#include <dust2/common.hpp>
#include <iostream>
#include <iterator>

// hardcoded as key to model structure
const int N_EPI_COMPARTMENTS = 4;
const int N_DATA_COMPARTMENTS = 1;
const int N_COMPARTMENTS = N_EPI_COMPARTMENTS + N_DATA_COMPARTMENTS;

// [[dust2::class(sirode)]]
// [[dust2::time_type(continuous)]]
// [[dust2::parameter(I0, constant = FALSE)]]
// [[dust2::parameter(N, constant = TRUE)]]
// [[dust2::parameter(beta, constant = FALSE)]]
// [[dust2::parameter(sigma, constant = FALSE)]]
// [[dust2::parameter(gamma, constant = FALSE)]]
// [[dust2::parameter(n_strata, constant = FALSE, type = "int")]]
class sirode {
public:
  sirode() = delete;

  using real_type = double;

  /// @brief Shared parameters and values.
  struct shared_state {
    real_type N;
    real_type I0;
    real_type beta;
    real_type sigma;
    real_type gamma;
  };

  /// @brief Internal state - unclear purpose.
  struct internal_state {
    double beta;
  };

  static internal_state build_internal(const shared_state &shared) {
    return internal_state{shared.beta};
  }

  /// @brief Holds incidence - unclear purpose.
  struct data_type {
    real_type incidence;
  };

  // unclear whether dust2/common.hpp links to monty - probably
  using rng_state_type = monty::random::generator<real_type>;

  /// @brief How compartments are packed.
  /// @param shared A `shared_state` object, unclear why needed.
  /// @return A custom packing specification object.
  static dust2::packing packing_state(const shared_state &shared) {
    // PASS STRATA DATA TO `shared` AND READ TO CREATE VECTORS FOR
    // PACKING VECTOR
    const std::vector<size_t> dim_vec(1, 1);
    const std::vector<size_t> dim_flag(1, 1);
    return dust2::packing{{"S", {}}, {"E", {}},         {"I", {}},
                          {"R", {}}, {"cases_inc", {}}, {"flag", {}}};
  }

  /// @brief Initialise shared parameters.
  /// @param pars A list of parameters passed from R.
  /// @return A shared parameters object.
  static shared_state build_shared(cpp11::list pars) {
    const real_type I0 = dust2::r::read_real(pars, "I0", 10);
    const real_type N = dust2::r::read_real(pars, "N", 1000);
    const real_type beta = dust2::r::read_real(pars, "beta", 0.2);
    const real_type sigma = dust2::r::read_real(pars, "sigma", 0.2);
    const real_type gamma = dust2::r::read_real(pars, "gamma", 0.1);
    return shared_state{N, I0, beta, sigma, gamma};
  }

  /// @brief Updated shared parameters.
  /// @param pars A list of parameters passed from R.
  /// @param shared A shared parameter object to update.
  static void update_shared(cpp11::list pars, shared_state &shared) {
    shared.I0 = dust2::r::read_real(pars, "I0", shared.I0);
    shared.beta = dust2::r::read_real(pars, "beta", shared.beta);
    shared.sigma = dust2::r::read_real(pars, "sigma", shared.sigma);
    shared.gamma = dust2::r::read_real(pars, "gamma", shared.gamma);
  }

  /// @brief Set initial values of the IVP model.
  /// @param time Time -- not used. Purpose unclear.
  /// @param shared Shared parameter object.
  /// @param internal Internal state object - not used. Purpose unclear.
  /// @param rng_state RNG state -- not used. Purpose unclear.
  /// @param state_next Next state as double value.
  static void initial(real_type time, const shared_state &shared,
                      internal_state &internal, rng_state_type &rng_state,
                      real_type *state_next) {
    // a response flag
    state_next[N_COMPARTMENTS] = 0.0;

    // map an Eigen container
    // TODO: figure out whether this is col or row major
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, N_COMPARTMENTS>> dx(
        &state_next[0], 1, N_COMPARTMENTS);

    // initially all zero, modify S and E
    dx.setZero();
    dx.col(0).setConstant(shared.N - shared.I0);
    dx.col(1).setConstant(shared.I0);
  }

  /// @brief RHS of the ODE model.
  /// @param time Time -- not used.
  /// @param state Pointer to state.
  /// @param shared Shared parameters.
  /// @param internal Internal state -- purpose unclear.
  /// @param state_deriv State change or dX.
  static void rhs(real_type time, const real_type *state,
                  const shared_state &shared, internal_state &internal,
                  real_type *state_deriv) {

    // map an Eigen container
    Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, N_COMPARTMENTS>> x(
        &state[0], 1, N_COMPARTMENTS);
    Eigen::Map<Eigen::Matrix<double, Eigen::Dynamic, N_COMPARTMENTS>> dx(
        &state_deriv[0], 1, N_COMPARTMENTS);
    // dx does not need to be set to zero as this is handled by zero_every()
    // seems like

    const auto rate_SE =
        shared.beta * x.col(0).array() * x.col(2).array() / shared.N;
    const auto rate_EI = shared.sigma * x.col(1).array();
    const auto rate_IR = shared.gamma * x.col(2).array();
    dx.col(0) = -rate_SE;
    dx.col(1) = rate_SE - rate_EI;
    dx.col(2) = rate_EI - rate_IR;
    dx.col(3) = rate_IR;
    dx.col(4) = rate_EI;
  }

  /// @brief Set every value to zero - unclear.
  /// @param shared Shared state -- unused.
  /// @return Probably an array of zeros.
  static auto zero_every(const shared_state &shared) {
    return dust2::zero_every_type<real_type>{
        {1, {4}}}; // zero only first three - works for single stratum
  }

  /// @brief Events for daedalus.
  /// @param shared Shared parameters.
  /// @param internal Intermediate containers.
  /// @return A container of events passed to the solver.
  static auto events(const shared_state &shared, internal_state &internal) {
    // check exposed and trigger event on root = 50 (prevalence)
    auto test = [&](double t, const double *y) {
      double diff = y[0] - 50.0; // check if root re: infected

      return diff;
    };

    auto action = [&](const double t, const double sign, double *y) {
      y[5] = 1.0; // set flag on
    };
    dust2::ode::event<real_type> e({1}, test, action,
                                   dust2::ode::root_type::increase);

    // check recovered and trigger event on root = 600 (epidemic size)
    auto test2 = [&](double t, const double *y) {
      double diff = y[0] - 600.0; // check if root re: infected

      return diff;
    };

    auto action2 = [&](const double t, const double sign, double *y) {
      y[5] = 0.0; // set flag off
    };
    dust2::ode::event<real_type> e2({3}, test, action);

    return dust2::ode::events_type<real_type>({e, e2});
  }
};
