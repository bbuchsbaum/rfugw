// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]

#include <RcppArmadillo.h>
#include "batch_worker.h"

// R adapter for the pure-C++ worker contract. No R object enters the native
// probe, and result wrapping occurs only after its worker region has joined.

// [[Rcpp::export]]
Rcpp::List cpp_thread_kernel_probe(
    int n_jobs = 16,
    int n_threads = 1,
    int nested_threads = 2,
    bool inject_failure = false) {
  const rfugw::BatchWorkerProbeResult out = rfugw::run_batch_worker_probe(
    n_jobs, n_threads, nested_threads, inject_failure
  );
  return Rcpp::List::create(
    Rcpp::Named("checksums") = out.checksums,
    Rcpp::Named("failed") = out.failed,
    Rcpp::Named("errors") = out.errors,
    Rcpp::Named("nested_suppressed") = out.nested_suppressed,
    Rcpp::Named("requested_threads") = out.requested_threads,
    Rcpp::Named("used_threads") = out.used_threads,
    Rcpp::Named("max_threads") = out.max_threads,
    Rcpp::Named("openmp") = out.max_threads > 1
  );
}
