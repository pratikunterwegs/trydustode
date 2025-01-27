#include <cpp11.hpp>
#include <cpp11eigen.hpp>
#include <dust2/common.hpp>
#include <iterator>

// [[dust2::class(sirode)]]
// [[dust2::time_type(continuous)]]
// [[dust2::has_compare()]]
// [[dust2::parameter(I0, constant = FALSE)]]
// [[dust2::parameter(N, constant = TRUE)]]
// [[dust2::parameter(beta, constant = FALSE)]]
// [[dust2::parameter(sigma, constant = FALSE)]]
// [[dust2::parameter(gamma, constant = FALSE)]]
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
    size_t n_strata;
  };

  /// @brief Internal state - unclear purpose.
  struct internal_state {};

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
    const std::vector<size_t> N_STRATA(1, shared.n_strata);
    return dust2::packing{{"S", N_STRATA},
                          {"E", N_STRATA},
                          {"I", N_STRATA},
                          {"R", N_STRATA},
                          {"cases_inc", N_STRATA}};
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
    const size_t n_strata = dust2::r::read_size(pars, "n_strata", 4);
    return shared_state{N, I0, beta, sigma, gamma, n_strata};
  }

  /// @brief Updated shared parameters.
  /// @param pars A list of parameters passed from R.
  /// @param shared A shared parameter object to update.
  static void update_shared(cpp11::list pars, shared_state &shared) {
    shared.I0 = dust2::r::read_real(pars, "I0", shared.I0);
    shared.beta = dust2::r::read_real(pars, "beta", shared.beta);
    shared.sigma = dust2::r::read_real(pars, "sigma", shared.sigma);
    shared.gamma = dust2::r::read_real(pars, "gamma", shared.gamma);
    shared.n_strata = dust2::r::read_real(pars, "n_strata", shared.n_strata);
  }

  /// @brief Return incidence data -- unclear purpose.
  /// @param r_data A list of R data to copy.
  /// @param shared Shared parameters -- unclear purpose.
  /// @return Data on incidence -- unclear purpose.
  static data_type build_data(cpp11::list r_data, const shared_state &shared) {
    auto data = static_cast<cpp11::list>(r_data);
    auto incidence = dust2::r::read_real(data, "incidence", NA_REAL);
    return data_type{incidence};
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
    // MAP AN EIGEN CONTAINER TO STATE_NEXT AND MODIFY CONTAINER VALUES
    state_next[0] = shared.N - shared.I0; // S
    state_next[1] = shared.N - shared.I0; // S
    state_next[2] = shared.I0;            // E
    state_next[3] = 0.0;                  // I
    state_next[4] = 0.0;                  // R
    state_next[5] = 0.0;                  // incidence
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
    // MAP EIGEN CONTAINER TO STATE AND D_STATE and read and modify
    // as in regular Eigen usage
    const auto S = state[0];
    const auto E = state[1];
    const auto I = state[2];
    const auto rate_SE = shared.beta * S * I / shared.N;
    const auto rate_EI = shared.sigma * E;
    const auto rate_IR = shared.gamma * I;
    state_deriv[0] = -rate_SE;
    state_deriv[1] = rate_SE;
    state_deriv[2] = rate_EI - rate_IR;
    state_deriv[3] = rate_IR;
    state_deriv[4] = rate_EI;

    // random calculation using cpp11eigen that goes nowhere
    // Eigen::Array<double, 10, 1> x;
    // x.setZero();

    // try mapping to state
    // NOTE: needs to be mapped to const and needs state
    Eigen::Map<const Eigen::VectorXd> y(state, 5);
    // y.setZero();
  }

  /// @brief Set every value to zero - unclear.
  /// @param shared Shared state -- unused.
  /// @return Probably an array of zeros.
  static auto zero_every(const shared_state &shared) {
    return dust2::zero_every_type<real_type>{{1, {3}}};
  }

  /// @brief Unclear what this does.
  /// @param time
  /// @param state
  /// @param data
  /// @param shared
  /// @param internal
  /// @param rng_state
  /// @return
  static real_type compare_data(const real_type time, const real_type *state,
                                const data_type &data,
                                const shared_state &shared,
                                internal_state &internal,
                                rng_state_type &rng_state) {
    const auto incidence_observed = data.incidence;
    if (std::isnan(incidence_observed)) {
      return 0;
    }
    const auto lambda = state[3];
    return monty::density::poisson(incidence_observed, lambda, true);
  }
};
