#include <dust2/common.hpp>

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
  static dust2::packing packing_state(const shared_state& shared) {
    return dust2::packing{{"S", {}}, {"E", {}}, {"I", {}}, {"R", {}}, {"cases_inc", {}}};
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
  static void update_shared(cpp11::list pars, shared_state& shared) {
    shared.I0 = dust2::r::read_real(pars, "I0", shared.I0);
    shared.beta = dust2::r::read_real(pars, "beta", shared.beta);
    shared.sigma = dust2::r::read_real(pars, "sigma", shared.sigma);
    shared.gamma = dust2::r::read_real(pars, "gamma", shared.gamma);
  }

  /// @brief Return incidence data -- unclear purpose.
  /// @param r_data A list of R data to copy.
  /// @param shared Shared parameters -- unclear purpose.
  /// @return Data on incidence -- unclear purpose.
  static data_type build_data(cpp11::list r_data, const shared_state& shared) {
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
  static void initial(real_type time,
                      const shared_state& shared,
                      internal_state& internal,
                      rng_state_type& rng_state,
                      real_type * state_next) {
    state_next[0] = shared.N - shared.I0;  // S
    state_next[1] = shared.I0;  // E
    state_next[2] = 0.0;  // I
    state_next[3] = 0.0;  // R
    state_next[4] = 0.0;  // incidence
  }

  /// @brief RHS of the ODE model.
  /// @param time Time -- not used.
  /// @param state Pointer to state.
  /// @param shared Shared parameters.
  /// @param internal Internal state -- purpose unclear.
  /// @param state_deriv State change or dX.
  static void rhs(real_type time,
                  const real_type * state,
                  const shared_state& shared,
                  internal_state& internal,
                  real_type * state_deriv) {
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
  }

  /// @brief Set every value to zero - unclear.
  /// @param shared Shared state -- unused.
  /// @return Probably an array of zeros.
  static auto zero_every(const shared_state& shared) {
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
  static real_type compare_data(const real_type time,
                                const real_type * state,
                                const data_type& data,
                                const shared_state& shared,
                                internal_state& internal,
                                rng_state_type& rng_state) {
    const auto incidence_observed = data.incidence;
    if (std::isnan(incidence_observed)) {
      return 0;
    }
    const auto lambda = state[3];
    return monty::density::poisson(incidence_observed, lambda, true);
  }
};
