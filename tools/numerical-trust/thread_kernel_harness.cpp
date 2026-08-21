#include "../../src/batch_worker.h"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

namespace {

bool same_checksums(
    const std::vector<double>& lhs,
    const std::vector<double>& rhs,
    double tolerance = 1e-12) {
  if (lhs.size() != rhs.size()) {
    return false;
  }
  for (std::size_t i = 0; i < lhs.size(); ++i) {
    const double scale = std::max(1.0, std::abs(lhs[i]));
    if (std::abs(lhs[i] - rhs[i]) > tolerance * scale) {
      return false;
    }
  }
  return true;
}

int failure_count(const std::vector<unsigned char>& failed) {
  int total = 0;
  for (unsigned char value : failed) {
    total += value == 0 ? 0 : 1;
  }
  return total;
}

}  // namespace

int main() {
  const auto serial = rfugw::run_batch_worker_probe(32, 1, 4, false);
  const auto two = rfugw::run_batch_worker_probe(32, 2, 4, false);
  const auto four = rfugw::run_batch_worker_probe(32, 4, 4, false);
  const auto injected = rfugw::run_batch_worker_probe(32, 4, 4, true);

  if (!same_checksums(serial.checksums, two.checksums) ||
      !same_checksums(serial.checksums, four.checksums)) {
    std::cerr << "serial/threaded checksum mismatch\n";
    return 1;
  }
  if (failure_count(injected.failed) != 1 ||
      injected.errors[16] != "injected worker failure") {
    std::cerr << "worker exception was not isolated and captured\n";
    return 2;
  }
  if (four.used_threads > 1) {
    for (unsigned char suppressed : four.nested_suppressed) {
      if (suppressed == 0) {
        std::cerr << "nested parallel work was not suppressed\n";
        return 3;
      }
    }
  }
  std::cout << "thread kernel probe passed: requested=4 used="
            << four.used_threads << " max=" << four.max_threads << "\n";
  return 0;
}
