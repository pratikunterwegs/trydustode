#include <dust2/common.hpp>

// [[dust2::class(bounce)]]
// [[dust2::time_type(continuous)]]
// [[dust2::parameter(height, rank = 0)]]
// [[dust2::parameter(velocity, rank = 0)]]
class bounce {
public:
  bounce() = delete;

  using real_type = double;

  struct shared_state {
    real_type g;
    real_type height;
    real_type velocity;
    real_type damp;
  };

  struct internal_state {};

  using rng_state_type = monty::random::generator<real_type>;

  static dust2::packing packing_state(const shared_state& shared) {
    return dust2::packing{{"height", {}}, {"velocity", {}}};
  }

  static void initial(real_type time,
                      const shared_state& shared,
                      internal_state& internal,
                      rng_state_type& rng_state,
                      real_type * state) {
    state[0] = shared.height;
    state[1] = shared.velocity;
  }

  static void rhs(real_type time,
                  const real_type * state,
                  const shared_state& shared,
                  internal_state& internal,
                  real_type * state_deriv) {
    state_deriv[0] = state[1];
    state_deriv[1] = -shared.g;
  }

  static shared_state build_shared(cpp11::list pars) {
    const real_type g = 9.8;
    const real_type height = dust2::r::read_real(pars, "height", 0);
    const real_type velocity = dust2::r::read_real(pars, "velocity", 10);
    const real_type damp = dust2::r::read_real(pars, "damp", 0.9);
    return shared_state{g, height, velocity, damp};
  }

  static void update_shared(cpp11::list pars, shared_state& shared) {
    shared.damp = dust2::r::read_real(pars, "damp", shared.damp);
  }

  static auto events(const shared_state& shared, internal_state& internal) {

    /// @brief Check if ball is at 5m.
    auto test = [&](const double t, const double *y) {
        double diff = 0.0;
        diff = y[0] - 5.0;

        return diff;
    };

    auto action = [&](const double t, const double sign, double* y) {
      y[0] = 5.0;
      y[1] = -shared.damp * y[1];
    };
    dust2::ode::event<real_type> e({0}, test, action);
    return dust2::ode::events_type<real_type>({e});
  }
};
