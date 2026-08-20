#pragma once  // certified exact-transport result contract

#include <armadillo>
#include <limits>
#include <stdexcept>
#include <string>

namespace rfugw {

enum class TransportTermination {
  optimal,
  max_iter,
  disconnected_basis,
  invalid_cycle,
  invalid_step,
  no_leaving_variable,
  numerical_failure
};

inline const char* transport_termination_name(TransportTermination reason) {
  switch (reason) {
  case TransportTermination::optimal:
    return "optimal";
  case TransportTermination::max_iter:
    return "max_iter";
  case TransportTermination::disconnected_basis:
    return "disconnected_basis";
  case TransportTermination::invalid_cycle:
    return "invalid_cycle";
  case TransportTermination::invalid_step:
    return "invalid_step";
  case TransportTermination::no_leaving_variable:
    return "no_leaving_variable";
  case TransportTermination::numerical_failure:
    return "numerical_failure";
  }
  return "numerical_failure";
}

enum class TransportTestFault {
  none,
  disconnected_basis,
  invalid_cycle,
  invalid_step,
  no_leaving_variable,
  numerical_failure
};

inline TransportTestFault parse_transport_test_fault(const std::string& fault) {
  if (fault.empty() || fault == "none") {
    return TransportTestFault::none;
  }
  if (fault == "disconnected_basis") {
    return TransportTestFault::disconnected_basis;
  }
  if (fault == "invalid_cycle") {
    return TransportTestFault::invalid_cycle;
  }
  if (fault == "invalid_step") {
    return TransportTestFault::invalid_step;
  }
  if (fault == "no_leaving_variable") {
    return TransportTestFault::no_leaving_variable;
  }
  if (fault == "numerical_failure") {
    return TransportTestFault::numerical_failure;
  }
  throw std::invalid_argument("unknown transport-simplex test fault: " + fault);
}

struct TransportSimplexResult {
  arma::mat plan;
  arma::vec source_potential;
  arma::vec target_potential;
  int iterations = 0;
  TransportTermination termination = TransportTermination::numerical_failure;
  bool certified = false;
  double primal_objective = std::numeric_limits<double>::quiet_NaN();
  double dual_objective = std::numeric_limits<double>::quiet_NaN();
  double duality_gap = std::numeric_limits<double>::quiet_NaN();
  double row_residual = std::numeric_limits<double>::infinity();
  double col_residual = std::numeric_limits<double>::infinity();
  double min_reduced_cost = -std::numeric_limits<double>::infinity();
  double feasibility_tolerance = std::numeric_limits<double>::quiet_NaN();
  double reduced_cost_tolerance = std::numeric_limits<double>::quiet_NaN();
  double duality_gap_tolerance = std::numeric_limits<double>::quiet_NaN();
  bool primal_feasible = false;
  bool dual_feasible = false;
  bool gap_certified = false;
};

}  // namespace rfugw
