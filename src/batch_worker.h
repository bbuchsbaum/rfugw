#ifndef RFUGW_BATCH_WORKER_H
#define RFUGW_BATCH_WORKER_H

#include <cstddef>
#include <algorithm>
#include <cmath>
#include <exception>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace rfugw {

// Execute one native batch item without permitting an exception to cross an
// OpenMP boundary. The caller owns one failure byte and one string per item, so
// simultaneous workers never write the same storage. R-facing conversion and
// error reporting deliberately remain outside this helper and outside workers.
template <typename Work>
inline void run_worker_guarded(
    std::size_t index,
    std::vector<unsigned char>& failed,
    std::vector<std::string>& errors,
    Work&& work) noexcept {
  try {
    std::forward<Work>(work)();
  } catch (const std::exception& e) {
    failed[index] = 1;
    errors[index] = e.what();
  } catch (...) {
    failed[index] = 1;
    errors[index] = "unknown native exception";
  }
}

struct BatchWorkerProbeResult {
  std::vector<double> checksums;
  std::vector<unsigned char> failed;
  std::vector<std::string> errors;
  std::vector<unsigned char> nested_suppressed;
  int requested_threads;
  int used_threads;
  int max_threads;
};

// Pure-C++ concurrency probe used both from R tests and a standalone TSan
// executable. Workers share one read-only cache, write only their indexed
// slots, suppress requested inner parallelism, and capture an injected error.
inline BatchWorkerProbeResult run_batch_worker_probe(
    int n_jobs,
    int n_threads,
    int nested_threads,
    bool inject_failure) {
  n_jobs = std::max(1, n_jobs);
  n_threads = std::max(1, n_threads);
  nested_threads = std::max(1, nested_threads);

  std::vector<double> cache(4096);
  for (std::size_t k = 0; k < cache.size(); ++k) {
    cache[k] = std::sin(static_cast<double>(k + 1) * 0.001);
  }

  BatchWorkerProbeResult result;
  result.checksums.assign(static_cast<std::size_t>(n_jobs), 0.0);
  result.failed.assign(static_cast<std::size_t>(n_jobs), 0);
  result.errors.resize(static_cast<std::size_t>(n_jobs));
  result.nested_suppressed.assign(static_cast<std::size_t>(n_jobs), 0);
  result.requested_threads = n_threads;
  result.used_threads = 1;
#ifdef _OPENMP
  result.max_threads = omp_get_max_threads();
  const int active_threads = std::min(n_threads, result.max_threads);
  std::vector<int> thread_ids(static_cast<std::size_t>(n_jobs), 0);
#pragma omp parallel for schedule(static) num_threads(active_threads) if(n_jobs > 1 && active_threads > 1)
  for (int i = 0; i < n_jobs; ++i) {
    const std::size_t idx = static_cast<std::size_t>(i);
    thread_ids[idx] = omp_get_thread_num();
    run_worker_guarded(idx, result.failed, result.errors, [&]() {
      if (nested_threads > 1 && omp_in_parallel() != 0) {
        result.nested_suppressed[idx] = 1;
      }
      if (inject_failure && i == n_jobs / 2) {
        throw std::runtime_error("injected worker failure");
      }
      double total = 0.0;
      for (std::size_t k = 0; k < cache.size(); ++k) {
        total += cache[k] * static_cast<double>((i + 1) * ((k % 17) + 1));
      }
      result.checksums[idx] = total;
    });
  }
  result.used_threads = 1 + *std::max_element(thread_ids.begin(), thread_ids.end());
#else
  result.max_threads = 1;
  for (int i = 0; i < n_jobs; ++i) {
    const std::size_t idx = static_cast<std::size_t>(i);
    run_worker_guarded(idx, result.failed, result.errors, [&]() {
      if (inject_failure && i == n_jobs / 2) {
        throw std::runtime_error("injected worker failure");
      }
      double total = 0.0;
      for (std::size_t k = 0; k < cache.size(); ++k) {
        total += cache[k] * static_cast<double>((i + 1) * ((k % 17) + 1));
      }
      result.checksums[idx] = total;
    });
  }
#endif
  return result;
}

}  // namespace rfugw

#endif
