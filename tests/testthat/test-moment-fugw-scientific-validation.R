moment_validation_test_env <- local({
  environment <- new.env(parent = asNamespace("rfugw"))
  sys.source(
    bench_test_resource("moment_fugw_scientific_validation.R"),
    envir = environment
  )
  environment
})

test_that("scientific validation protocol was frozen before evaluation", {
  protocol <- moment_validation_test_env$moment_validation_protocol()

  expect_equal(protocol$schema_version, 2)
  expect_identical(protocol$status, "frozen_before_evaluation")
  expect_false(protocol$pre_evaluation_amendment$
               evaluation_seed_accessed_before_freeze)
  expect_false(protocol$post_evaluation_replay_hardening$
               estimand_threshold_or_evaluation_result_changed)
  expect_true(protocol$post_evaluation_replay_hardening$
              semantic_fixture_fingerprint_preserved)
  expect_length(
    intersect(protocol$image$tuning_seeds, protocol$image$evaluation_seeds),
    0L
  )
  expect_length(
    intersect(protocol$surface$tuning_seeds, protocol$surface$evaluation_seeds),
    0L
  )
  expect_equal(protocol$image$overlap_percent, c(100, 80, 60, 40))
  expect_equal(
    protocol$image$mask_kind,
    c("none", "crop", "lesion", "crop_plus_lesion")
  )
  expect_identical(
    protocol$failure_policy,
    "retain_every_requested_case_with_solver_success_and_failure_reason"
  )
  expect_identical(protocol$surface$dataset, "neuroatlas::fsaverage6")
  expect_equal(protocol$surface$landmark_count, 36)
  expect_equal(protocol$surface$dense_parity_landmark_count, 16)
  expect_identical(
    protocol$surface$dense_parity_fixture_fingerprint_md5,
    "ebaaa6eecd52a7dd7ee339fe8d294cbd"
  )
  expect_equal(protocol$surface$train_anchor_ids, c(1, 11, 21, 31))
  expect_equal(protocol$surface$heldout_anchor_ids, c(6, 16, 26))
})

test_that("scientific provenance treats unavailable Git state as unknown", {
  outside <- tempfile("moment-science-no-git-")
  dir.create(outside)
  provenance <- moment_validation_test_env$moment_git_provenance(outside)

  expect_true(is.na(provenance$git_commit))
  expect_true(is.na(provenance$git_dirty))
  expect_true(is.na(provenance$git_status_entry_count))
  expect_false(provenance$git_provenance_complete)
  expect_false(identical(provenance$git_dirty, FALSE))
})

test_that("planted image fixtures encode deformation, overlap, and heldout maps", {
  protocol <- moment_validation_test_env$moment_validation_protocol()$image
  first <- moment_validation_test_env$moment_make_image_case(4101L, 60, protocol)
  second <- moment_validation_test_env$moment_make_image_case(4101L, 60, protocol)

  expect_identical(first, second)
  expect_equal(nrow(first$source_coordinates), 100L)
  expect_equal(nrow(first$target_coordinates), 60L)
  expect_length(first$overlap, 60L)
  expect_length(first$occluded, 40L)
  expect_identical(first$mask_kind, "lesion")
  expect_equal(ncol(first$source_features), 4L)
  expect_equal(ncol(first$source_heldout), 3L)
  expect_gt(stats::sd(first$source_weights), 0)
  expect_gt(stats::sd(first$target_weights), 0)
  undeformed <- first$source_coordinates[first$target_order, , drop = FALSE]
  expect_gt(mean(abs(first$target_coordinates - undeformed)), 0.01)

  severe <- moment_validation_test_env$moment_make_image_case(
    4101L, 40, protocol
  )
  expect_identical(severe$mask_kind, "crop_plus_lesion")
  expect_length(severe$occluded, 60L)
})

test_that("retained image evaluation is complete, paired, and passes margins", {
  protocol <- moment_validation_test_env$moment_validation_protocol()$image
  evidence <- utils::read.csv(
    bench_test_resource("moment-fugw-image-evidence.csv"),
    stringsAsFactors = FALSE
  )

  expect_equal(nrow(evidence), 40L)
  expect_setequal(unique(evidence$seed), protocol$evaluation_seeds)
  expect_setequal(unique(evidence$overlap_percent), protocol$overlap_percent)
  expect_setequal(unique(evidence$method), c("unbalanced_fugw", "balanced_fgw"))
  expect_true(all(table(evidence$method, evidence$overlap_percent) == 5L))
  expect_true(all(evidence$solver_success))
  expect_true(all(evidence$evaluation_split == "frozen_evaluation"))
  expect_true(all(evidence$suite_gate))
  expect_true(all(evidence$protocol_schema_version == 2L))
  expect_true(all(evidence$protocol_status == "frozen_before_evaluation"))
  expect_true(all(!evidence$validation_materialization[
    evidence$method == "unbalanced_fugw"
  ]))
  expect_true(all(evidence$validation_materialization[
    evidence$method == "balanced_fgw"
  ]))

  metrics <- c(
    "transported_mass", "recovered_overlap_mass", "discarded_source_mass",
    "discarded_target_mass", "false_occluded_source_mass",
    "weighted_endpoint_error", "heldout_rmse", "heldout_correlation",
    "local_distortion"
  )
  expect_true(all(vapply(evidence[metrics], function(x) all(is.finite(x)),
                         logical(1))))
  summary <- moment_validation_test_env$moment_image_summary(evidence, protocol)
  expect_equal(summary$overlap_percent, c(80, 60, 40))
  expect_true(all(summary$paired_complete))
  expect_true(all(summary$paired_success_count == 5L))
  expect_true(all(summary$false_mass_gate))
  expect_true(all(summary$endpoint_gate))
  expect_true(all(summary$recovered_mass_gate))
  expect_true(all(
    summary$unbalanced_median_false_mass < summary$balanced_median_false_mass
  ))
  expect_lte(
    max(summary$endpoint_error_ratio),
    protocol$partial_overlap_gate$median_endpoint_error_ratio_max
  )
})

test_that("retained cortical evaluation and geometry curve are certified", {
  protocol <- moment_validation_test_env$moment_validation_protocol()$surface
  evidence <- utils::read.csv(
    bench_test_resource("moment-fugw-surface-evidence.csv"),
    stringsAsFactors = FALSE
  )
  curve <- utils::read.csv(
    bench_test_resource("moment-fugw-geometry-rank-evidence.csv"),
    stringsAsFactors = FALSE
  )

  expect_equal(nrow(evidence), 5L)
  expect_setequal(unique(evidence$seed), protocol$evaluation_seeds)
  expect_true(all(evidence$solver_success))
  expect_true(all(evidence$geometry_certified))
  expect_true(all(evidence$geometry_exact))
  expect_true(all(evidence$geometry_relative_error == 0))
  expect_true(all(evidence$geometry_heldout_error == 0))
  expect_true(all(evidence$dataset == protocol$dataset))
  expect_true(all(evidence$dataset_package_version ==
                  protocol$dataset_package_version))
  expect_true(all(evidence$fixture_fingerprint_md5 ==
                  protocol$fixture_fingerprint_md5))
  expect_true(all(evidence$full_vertex_count == protocol$full_vertex_count))
  expect_true(all(evidence$full_face_count == protocol$full_face_count))
  expect_true(all(evidence$n_source == protocol$landmark_count))
  expect_true(all(evidence$area_weight_method == protocol$area_weight_method))
  expect_true(all(evidence$functional_map_source ==
                  protocol$functional_map_source))
  expect_true(all(evidence$evaluation_split == "frozen_evaluation"))
  expect_true(all(evidence$suite_gate))
  expect_true(all(evidence$dense_parity_passed))
  expect_lte(
    max(evidence$dense_parity_maximum_error),
    protocol$dense_parity_tolerance
  )
  surface_metrics <- c(
    "heldout_correlation", "heldout_r2",
    "expected_geodesic_endpoint_error", "transported_mass",
    "discarded_source_mass", "discarded_target_mass",
    "coupling_l1_discrepancy", "coupling_max_discrepancy"
  )
  expect_true(all(vapply(evidence[surface_metrics],
                         function(x) all(is.finite(x)), logical(1))))

  expect_equal(nrow(curve), length(protocol$embedding_curve$ranks))
  expect_equal(curve$embedding_rank, protocol$embedding_curve$ranks)
  expect_true(all(curve$solver_success))
  expect_true(all(curve$geometry_certified))
  expect_true(all(!curve$geometry_exact))
  expect_true(all(is.finite(curve$geometry_relative_error)))
  expect_true(all(is.finite(curve$geometry_heldout_error)))
  expect_true(all(is.finite(curve$heldout_correlation)))
  expect_true(all(curve$numerical_and_geometry_error_separate))
  expect_gt(min(curve$geometry_relative_error),
            max(curve$coupling_l1_discrepancy))
})

test_that("pinned cortical fixture and dense parity can be replayed", {
  protocol <- moment_validation_test_env$moment_validation_protocol()$surface

  first <- moment_validation_test_env$moment_make_surface_case(
    5101L, protocol$landmark_count
  )
  second <- moment_validation_test_env$moment_make_surface_case(
    5101L, protocol$landmark_count
  )
  expect_identical(first, second)
  expect_invisible(
    moment_validation_test_env$.moment_validate_surface_fixture(first, protocol)
  )
  expect_equal(sum(first$source_weights), 1, tolerance = 1e-14)
  expect_gt(stats::sd(first$source_weights), 0)
  expect_equal(dim(first$source_geodesic), c(36L, 36L))
  expect_true(all(is.finite(first$source_geodesic)))
  expect_equal(ncol(first$source_features), 4L)
  expect_equal(ncol(first$source_heldout), 3L)

  tampered <- protocol
  tampered$fixture_fingerprint_md5 <- "tampered"
  expect_error(
    moment_validation_test_env$.moment_validate_surface_fixture(first, tampered),
    "fixture_fingerprint"
  )
  parity <- moment_validation_test_env$moment_surface_dense_parity(protocol)
  expect_true(parity$factorized_certified)
  expect_true(parity$dense_converged)
  expect_true(parity$fixture_fingerprint_passed)
  expect_identical(
    parity$fixture_fingerprint_md5,
    protocol$dense_parity_fixture_fingerprint_md5
  )
  expect_true(parity$passed)
  expect_lte(parity$maximum_parity_error, protocol$dense_parity_tolerance)
})
