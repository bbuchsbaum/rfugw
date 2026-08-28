moment_validation_resource <- function(name) {
  installed <- system.file("bench", name, package = "rfugw")
  if (nzchar(installed) && file.exists(installed)) return(installed)
  source_path <- file.path("inst", "bench", name)
  if (file.exists(source_path)) return(source_path)
  stop("Moment-FUGW validation resource not found: ", name, call. = FALSE)
}

moment_validation_protocol <- function(path = NULL) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Scientific validation requires jsonlite.", call. = FALSE)
  }
  if (is.null(path)) {
    path <- moment_validation_resource("moment-fugw-validation-protocol.json")
  }
  jsonlite::fromJSON(path, simplifyVector = TRUE)
}

moment_git_provenance <- function(repository_root = ".") {
  repository_root <- tryCatch(
    normalizePath(repository_root, mustWork = TRUE),
    error = function(e) NA_character_
  )
  capture <- function(arguments) {
    if (is.na(repository_root)) {
      return(list(ok = FALSE, status = NA_integer_, value = character()))
    }
    value <- suppressWarnings(tryCatch(
      system2(
        "git", c("-C", repository_root, arguments),
        stdout = TRUE, stderr = FALSE
      ),
      error = function(e) structure(character(), status = 127L)
    ))
    status <- attr(value, "status", exact = TRUE)
    if (is.null(status)) status <- 0L
    list(
      ok = identical(as.integer(status), 0L),
      status = as.integer(status),
      value = unname(as.character(value))
    )
  }
  commit_result <- capture(c("rev-parse", "--verify", "HEAD"))
  status_result <- capture(c(
    "status", "--porcelain", "--untracked-files=normal"
  ))
  commit <- if (isTRUE(commit_result$ok) &&
      length(commit_result$value) == 1L &&
      grepl("^[[:xdigit:]]{40}$", commit_result$value[[1L]])) {
    tolower(commit_result$value[[1L]])
  } else {
    NA_character_
  }
  status_ok <- isTRUE(status_result$ok)
  list(
    repository_root = repository_root,
    git_commit = commit,
    git_dirty = if (status_ok) length(status_result$value) > 0L else NA,
    git_status_entry_count = if (status_ok) {
      as.integer(length(status_result$value))
    } else {
      NA_integer_
    },
    git_provenance_complete = !is.na(commit) && status_ok,
    git_commit_command_status = commit_result$status,
    git_status_command_status = status_result$status
  )
}

.moment_weighted_mean <- function(x, weights) {
  sum(weights * x) / sum(weights)
}

.moment_weighted_rmse <- function(value, truth, weights) {
  difference <- rowMeans((value - truth)^2)
  sqrt(.moment_weighted_mean(difference, weights))
}

.moment_weighted_correlation <- function(value, truth, weights) {
  x <- as.numeric(value)
  y <- as.numeric(truth)
  w <- rep(weights, ncol(value))
  w <- w / sum(w)
  x_mean <- sum(w * x)
  y_mean <- sum(w * y)
  covariance <- sum(w * (x - x_mean) * (y - y_mean))
  denominator <- sqrt(
    sum(w * (x - x_mean)^2) * sum(w * (y - y_mean)^2)
  )
  if (denominator <= .Machine$double.eps) return(NA_real_)
  covariance / denominator
}

.moment_weighted_r2 <- function(value, truth, weights) {
  w <- rep(weights, ncol(value))
  y <- as.numeric(truth)
  x <- as.numeric(value)
  center <- sum(w * y) / sum(w)
  1 - sum(w * (y - x)^2) / max(
    sum(w * (y - center)^2), .Machine$double.eps
  )
}

.moment_normalized_apply <- function(plan, target_values) {
  source_mass <- transport_plan_mass(plan, "source")
  prediction <- transport_plan_apply(plan, target_values)
  active <- source_mass > sqrt(.Machine$double.eps)
  prediction[active, ] <- prediction[active, , drop = FALSE] /
    source_mass[active]
  list(prediction = prediction, source_mass = source_mass, active = active)
}

.moment_image_latent_channels <- function(coordinates) {
  x <- coordinates[, 1L]
  y <- coordinates[, 2L]
  train <- cbind(
    exp(-((x - .25)^2 + (y - .3)^2) / .045),
    exp(-((x - .72)^2 + (y - .68)^2) / .06),
    sin(pi * x) * cos(pi * y),
    x - .6 * y
  )
  heldout <- cbind(
    cos(2 * pi * x) + .25 * y,
    sin(2 * pi * y) - .2 * x,
    exp(-((x - .5)^2 + (y - .5)^2) / .08)
  )
  list(train = train, heldout = heldout)
}

.moment_image_overlap <- function(coordinates, overlap_percent, mask_kind) {
  n <- nrow(coordinates)
  keep_n <- as.integer(round(n * overlap_percent / 100))
  x <- coordinates[, 1L]
  y <- coordinates[, 2L]
  if (identical(mask_kind, "none")) return(seq_len(n))
  if (identical(mask_kind, "crop")) {
    return(order(x, y)[seq_len(keep_n)])
  }
  radius <- (x - .5)^2 + (y - .5)^2
  if (identical(mask_kind, "lesion")) {
    return(order(radius, decreasing = TRUE)[seq_len(keep_n)])
  }
  crop <- order(x, y)[seq_len(min(n, keep_n + 10L))]
  crop[order(radius[crop], decreasing = TRUE)[seq_len(keep_n)]]
}

moment_make_image_case <- function(seed, overlap_percent, protocol = NULL) {
  if (is.null(protocol)) protocol <- moment_validation_protocol()$image
  set.seed(seed + as.integer(overlap_percent) * 1000L)
  side <- as.integer(protocol$grid_side)
  axis <- seq(0, 1, length.out = side)
  source_coordinates <- as.matrix(expand.grid(x = axis, y = axis))
  mask_position <- match(overlap_percent, protocol$overlap_percent)
  mask_kind <- protocol$mask_kind[[mask_position]]
  overlap <- sort(.moment_image_overlap(
    source_coordinates, overlap_percent, mask_kind
  ))
  occluded <- setdiff(seq_len(nrow(source_coordinates)), overlap)
  target_order <- sample(overlap)
  latent <- .moment_image_latent_channels(source_coordinates)
  noise_sd <- protocol$feature_noise_sd
  source_features <- latent$train + matrix(
    rnorm(length(latent$train), sd = noise_sd),
    nrow(latent$train), ncol(latent$train)
  )
  target_features <- latent$train[target_order, , drop = FALSE] + matrix(
    rnorm(length(target_order) * ncol(latent$train), sd = noise_sd),
    length(target_order), ncol(latent$train)
  )
  target_heldout <- latent$heldout[target_order, , drop = FALSE] + matrix(
    rnorm(length(target_order) * ncol(latent$heldout), sd = noise_sd / 4),
    length(target_order), ncol(latent$heldout)
  )
  target_coordinates <- source_coordinates[target_order, , drop = FALSE]
  x <- target_coordinates[, 1L]
  y <- target_coordinates[, 2L]
  target_coordinates <- cbind(
    x + .08 * sin(pi * y),
    y + .06 * x * (1 - x)
  )
  angle <- .12
  rotation <- matrix(
    c(cos(angle), -sin(angle), sin(angle), cos(angle)),
    2L, 2L, byrow = TRUE
  )
  target_coordinates <- target_coordinates %*% t(rotation)
  raw_source_weights <- .55 + .7 * source_coordinates[, 1L] +
    .15 * cos(2 * pi * source_coordinates[, 2L])
  source_weights <- raw_source_weights / sum(raw_source_weights)
  raw_target_weights <- raw_source_weights[target_order] *
    (1 + .12 * target_coordinates[, 2L])
  target_weights <- raw_target_weights / sum(raw_target_weights)
  truth_target_index <- match(seq_len(nrow(source_coordinates)), target_order)
  list(
    seed = seed,
    overlap_percent = overlap_percent,
    mask_kind = mask_kind,
    source_coordinates = source_coordinates,
    target_coordinates = target_coordinates,
    source_features = source_features,
    target_features = target_features,
    source_heldout = latent$heldout,
    target_heldout = target_heldout,
    source_weights = source_weights,
    target_weights = target_weights,
    overlap = overlap,
    occluded = occluded,
    target_order = target_order,
    truth_target_index = truth_target_index,
    feature_noise_sd = noise_sd
  )
}

.moment_image_metrics <- function(
    plan, problem, method, solver_success, status, failure_reason = "") {
  empty <- data.frame(
    method = method,
    solver_success = solver_success,
    status = status,
    failure_reason = failure_reason,
    transported_mass = NA_real_,
    recovered_overlap_mass = NA_real_,
    recovered_overlap_fraction = NA_real_,
    discarded_source_mass = NA_real_,
    discarded_target_mass = NA_real_,
    false_occluded_source_mass = NA_real_,
    false_occluded_mass_fraction = NA_real_,
    weighted_endpoint_error = NA_real_,
    heldout_rmse = NA_real_,
    heldout_correlation = NA_real_,
    local_distortion = NA_real_,
    stringsAsFactors = FALSE
  )
  if (!solver_success) return(empty)
  applied_coordinates <- .moment_normalized_apply(
    plan, problem$target_coordinates
  )
  applied_heldout <- .moment_normalized_apply(plan, problem$target_heldout)
  source_mass <- applied_coordinates$source_mass
  target_mass <- transport_plan_mass(plan, "target")
  transported_mass <- sum(source_mass)
  overlap <- problem$overlap
  truth_coordinates <- problem$target_coordinates[
    problem$truth_target_index[overlap], , drop = FALSE
  ]
  endpoint <- sqrt(rowSums(
    (applied_coordinates$prediction[overlap, , drop = FALSE] -
       truth_coordinates)^2
  ))
  active_overlap <- overlap[
    applied_coordinates$active[overlap] & applied_heldout$active[overlap]
  ]
  overlap_weights <- problem$source_weights[active_overlap]
  predicted_coordinates <- applied_coordinates$prediction[
    active_overlap, , drop = FALSE
  ]
  true_coordinates <- problem$target_coordinates[
    problem$truth_target_index[active_overlap], , drop = FALSE
  ]
  predicted_distances <- as.matrix(stats::dist(predicted_coordinates))
  true_distances <- as.matrix(stats::dist(true_coordinates))
  pair_weights <- outer(overlap_weights, overlap_weights)
  local_distortion <- sqrt(
    sum(pair_weights * (predicted_distances - true_distances)^2) /
      max(sum(pair_weights * true_distances^2), .Machine$double.eps)
  )
  recovered_overlap_mass <- sum(source_mass[overlap])
  false_mass <- sum(source_mass[problem$occluded])
  data.frame(
    method = method,
    solver_success = TRUE,
    status = status,
    failure_reason = "",
    transported_mass = transported_mass,
    recovered_overlap_mass = recovered_overlap_mass,
    recovered_overlap_fraction = recovered_overlap_mass /
      sum(problem$source_weights[overlap]),
    discarded_source_mass = max(0, sum(problem$source_weights) - transported_mass),
    discarded_target_mass = max(0, sum(problem$target_weights) - transported_mass),
    false_occluded_source_mass = false_mass,
    false_occluded_mass_fraction = false_mass / max(
      transported_mass, .Machine$double.eps
    ),
    weighted_endpoint_error = .moment_weighted_mean(
      endpoint, problem$source_weights[overlap]
    ),
    heldout_rmse = .moment_weighted_rmse(
      applied_heldout$prediction[active_overlap, , drop = FALSE],
      problem$source_heldout[active_overlap, , drop = FALSE],
      overlap_weights
    ),
    heldout_correlation = .moment_weighted_correlation(
      applied_heldout$prediction[active_overlap, , drop = FALSE],
      problem$source_heldout[active_overlap, , drop = FALSE],
      overlap_weights
    ),
    local_distortion = local_distortion,
    stringsAsFactors = FALSE
  )
}

moment_run_image_case <- function(seed, overlap_percent, protocol = NULL) {
  if (is.null(protocol)) protocol <- moment_validation_protocol()$image
  problem <- moment_make_image_case(seed, overlap_percent, protocol)
  Cx <- sqeuclidean_cost(problem$source_coordinates)
  Cy <- sqeuclidean_cost(problem$target_coordinates)
  M <- sqeuclidean_cost(problem$source_features, problem$target_features)
  unbalanced <- tryCatch({
    control <- protocol$unbalanced
    fit <- fugw_factorized(
      Cx, Cy,
      wx = problem$source_weights,
      wy = problem$target_weights,
      M = M,
      feature_weight = control$feature_weight,
      structure_weight = control$structure_weight,
      reg_marginals = control$reg_marginals,
      epsilon = control$epsilon,
      max_iter = as.integer(control$max_iter),
      tol = control$tol,
      max_iter_ot = as.integer(control$max_iter_ot),
      tol_ot = control$tol_ot,
      block_size = as.integer(control$block_size)
    )
    success <- isTRUE(fit$final_stationarity_certified) &&
      isTRUE(fit$geometry_certified)
    metric <- .moment_image_metrics(
      fit$plans$sample, problem, "unbalanced_fugw", success, fit$status,
      if (success) "" else paste(fit$certification_failures, collapse = ";")
    )
    metric$coupling_l1_discrepancy <- fit$coupling_discrepancy$l1
    metric$geometry_relative_error <- fit$geometry_relative_error
    metric
  }, error = function(error) {
    metric <- .moment_image_metrics(
      NULL, problem, "unbalanced_fugw", FALSE, "error", conditionMessage(error)
    )
    metric$coupling_l1_discrepancy <- NA_real_
    metric$geometry_relative_error <- NA_real_
    metric
  })
  balanced <- tryCatch({
    control <- protocol$balanced_comparator
    fit <- fgw_entropic(
      M = cost_block(M),
      C1 = cost_block(Cx),
      C2 = cost_block(Cy),
      p = problem$source_weights,
      q = problem$target_weights,
      feature_weight = control$feature_weight,
      structure_weight = control$structure_weight,
      epsilon = control$epsilon,
      max_iter = as.integer(control$max_iter),
      tol = control$tol,
      sinkhorn_max_iter = as.integer(control$sinkhorn_max_iter),
      sinkhorn_tol = control$sinkhorn_tol,
      precision = control$precision,
      sinkhorn_method = control$sinkhorn_method
    )
    success <- isTRUE(fit$converged) && isTRUE(fit$inner_converged) &&
      isTRUE(fit$feasible)
    metric <- .moment_image_metrics(
      as_transport_plan(fit$plan), problem, "balanced_fgw", success,
      fit$status, if (success) "" else fit$termination_reason
    )
    metric$coupling_l1_discrepancy <- 0
    metric$geometry_relative_error <- 0
    metric
  }, error = function(error) {
    metric <- .moment_image_metrics(
      NULL, problem, "balanced_fgw", FALSE, "error", conditionMessage(error)
    )
    metric$coupling_l1_discrepancy <- NA_real_
    metric$geometry_relative_error <- 0
    metric
  })
  rows <- rbind(unbalanced, balanced)
  rows$seed <- seed
  rows$overlap_percent <- overlap_percent
  rows$mask_kind <- problem$mask_kind
  rows$feature_noise_sd <- problem$feature_noise_sd
  rows$n_source <- nrow(problem$source_coordinates)
  rows$n_target <- nrow(problem$target_coordinates)
  rows$unequal_weights <- TRUE
  rows$validation_materialization <- rows$method == "balanced_fgw"
  rows
}

moment_image_summary <- function(rows, protocol = NULL) {
  if (is.null(protocol)) protocol <- moment_validation_protocol()$image
  requested_overlaps <- sort(
    protocol$overlap_percent[protocol$overlap_percent < 100],
    decreasing = TRUE
  )
  summaries <- do.call(rbind, lapply(
    requested_overlaps,
    function(overlap) {
      current <- rows[rows$overlap_percent == overlap, , drop = FALSE]
      u <- current[current$method == "unbalanced_fugw", , drop = FALSE]
      b <- current[current$method == "balanced_fgw", , drop = FALSE]
      paired_seeds <- intersect(
        u$seed[u$solver_success], b$seed[b$solver_success]
      )
      requested_seeds <- union(u$seed, b$seed)
      u <- u[u$seed %in% paired_seeds, , drop = FALSE]
      b <- b[b$seed %in% paired_seeds, , drop = FALSE]
      paired_complete <- length(paired_seeds) == length(requested_seeds) &&
        nrow(u) == length(requested_seeds) &&
        nrow(b) == length(requested_seeds)
      data.frame(
        overlap_percent = overlap,
        requested_case_count = length(requested_seeds),
        paired_success_count = length(paired_seeds),
        paired_complete = paired_complete,
        unbalanced_median_false_mass = if (nrow(u)) {
          median(u$false_occluded_source_mass)
        } else NA_real_,
        balanced_median_false_mass = if (nrow(b)) {
          median(b$false_occluded_source_mass)
        } else NA_real_,
        endpoint_error_ratio = if (nrow(u) && nrow(b)) {
          median(u$weighted_endpoint_error) /
            median(b$weighted_endpoint_error)
        } else NA_real_,
        unbalanced_median_recovered_overlap_fraction =
          if (nrow(u)) median(u$recovered_overlap_fraction) else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  ))
  gate <- protocol$partial_overlap_gate
  summaries$false_mass_gate <- summaries$paired_complete &
    is.finite(summaries$unbalanced_median_false_mass) &
    is.finite(summaries$balanced_median_false_mass) &
    summaries$unbalanced_median_false_mass <
      summaries$balanced_median_false_mass
  summaries$endpoint_gate <- summaries$paired_complete &
    is.finite(summaries$endpoint_error_ratio) &
    summaries$endpoint_error_ratio <= gate$median_endpoint_error_ratio_max
  summaries$recovered_mass_gate <-
    summaries$paired_complete &
      is.finite(summaries$unbalanced_median_recovered_overlap_fraction) &
      summaries$unbalanced_median_recovered_overlap_fraction >=
      gate$median_recovered_overlap_fraction_min
  summaries
}

moment_run_image_suite <- function(
    seeds = NULL,
    output = file.path("inst", "bench", "moment-fugw-image-evidence.csv")) {
  protocol <- moment_validation_protocol()$image
  if (is.null(seeds)) seeds <- protocol$evaluation_seeds
  rows <- do.call(rbind, lapply(seeds, function(seed) {
    do.call(rbind, lapply(protocol$overlap_percent, function(overlap) {
      message("Image validation seed ", seed, ", overlap ", overlap)
      moment_run_image_case(seed, overlap, protocol)
    }))
  }))
  summary <- moment_image_summary(rows, protocol)
  rows$evaluation_split <- ifelse(
    rows$seed %in% protocol$evaluation_seeds, "frozen_evaluation", "tuning"
  )
  summary_index <- match(rows$overlap_percent, summary$overlap_percent)
  for (column in setdiff(names(summary), "overlap_percent")) {
    rows[[column]] <- summary[[column]][summary_index]
  }
  rows$suite_gate <- all(rows$solver_success) &&
    all(summary$false_mass_gate) && all(summary$endpoint_gate) &&
    all(summary$recovered_mass_gate)
  full_protocol <- moment_validation_protocol()
  rows$protocol_schema_version <- full_protocol$schema_version
  rows$protocol_status <- full_protocol$status
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(rows, output, row.names = FALSE)
  attr(rows, "summary") <- summary
  rows
}

.moment_validation_cache <- new.env(parent = emptyenv())

.moment_weighted_tabulate <- function(group, weights, nbins) {
  sums <- rowsum(weights, group = group, reorder = FALSE)
  result <- numeric(as.integer(nbins))
  result[as.integer(rownames(sums))] <- sums[, 1L]
  result
}

.moment_triangle_areas <- function(vertices, faces) {
  ab <- vertices[faces[, 2L], , drop = FALSE] -
    vertices[faces[, 1L], , drop = FALSE]
  ac <- vertices[faces[, 3L], , drop = FALSE] -
    vertices[faces[, 1L], , drop = FALSE]
  normal <- cbind(
    ab[, 2L] * ac[, 3L] - ab[, 3L] * ac[, 2L],
    ab[, 3L] * ac[, 1L] - ab[, 1L] * ac[, 3L],
    ab[, 1L] * ac[, 2L] - ab[, 2L] * ac[, 1L]
  )
  sqrt(rowSums(normal^2)) / 2
}

.moment_farthest_landmarks <- function(vertices, count) {
  count <- as.integer(count)
  if (count < 8L || count > nrow(vertices)) {
    stop("landmark_count must be between 8 and the full mesh size.", call. = FALSE)
  }
  selected <- integer(count)
  selected[[1L]] <- which.max(vertices[, 3L])
  nearest_squared_distance <- rep(Inf, nrow(vertices))
  for (index in seq_len(count)) {
    landmark <- vertices[selected[[index]], ]
    distance <- rowSums(sweep(vertices, 2L, landmark, "-")^2)
    nearest_squared_distance <- pmin(nearest_squared_distance, distance)
    nearest_squared_distance[selected[seq_len(index)]] <- -Inf
    if (index < count) {
      selected[[index + 1L]] <- which.max(nearest_squared_distance)
    }
  }
  selected
}

.moment_aggregate_area_weights <- function(
    vertices, landmarks, vertex_areas, block_size = 2048L) {
  landmark_vertices <- vertices[landmarks, , drop = FALSE]
  landmark_norms <- rowSums(landmark_vertices^2)
  assignment <- integer(nrow(vertices))
  starts <- seq.int(1L, nrow(vertices), by = as.integer(block_size))
  for (start in starts) {
    rows <- start:min(nrow(vertices), start + block_size - 1L)
    current <- vertices[rows, , drop = FALSE]
    squared_distance <- outer(rowSums(current^2), landmark_norms, "+") -
      2 * tcrossprod(current, landmark_vertices)
    assignment[rows] <- max.col(-squared_distance, ties.method = "first")
  }
  weights <- .moment_weighted_tabulate(
    assignment, vertex_areas, nbins = length(landmarks)
  )
  weights / sum(weights)
}

.moment_fixture_fingerprint <- function(fixture) {
  numeric_text <- function(x) {
    format(
      as.numeric(x), digits = 17L, scientific = TRUE, trim = TRUE,
      decimal.mark = "."
    )
  }
  payload <- c(
    "moment_fsaverage6_lh_pial_fixture_v1",
    fixture$dataset,
    fixture$hemisphere,
    fixture$surface,
    fixture$landmark_method,
    fixture$area_weight_method,
    paste(fixture$landmark_vertex_ids, collapse = ","),
    numeric_text(t(fixture$vertices)),
    numeric_text(fixture$weights),
    numeric_text(t(fixture$geodesic))
  )
  path <- tempfile("moment-fugw-fixture-")
  connection <- file(path, open = "wb")
  connection_open <- TRUE
  on.exit({
    if (connection_open) close(connection)
    unlink(path)
  }, add = TRUE)
  writeBin(charToRaw(paste(payload, collapse = "\n")), connection)
  close(connection)
  connection_open <- FALSE
  unname(tools::md5sum(path))
}

.moment_cortical_fixture <- function(landmark_count = 36L) {
  landmark_count <- as.integer(landmark_count)
  cache_key <- paste0("fsaverage6_lh_pial_retained_", landmark_count)
  if (exists(cache_key, envir = .moment_validation_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .moment_validation_cache, inherits = FALSE))
  }
  if (!landmark_count %in% c(16L, 36L)) {
    stop(
      "Only the retained 16- and 36-landmark cortical fixtures are available.",
      call. = FALSE
    )
  }
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("Cortical fixture validation requires jsonlite.", call. = FALSE)
  }
  provenance_path <- moment_validation_resource(
    "moment-fugw-fsaverage6-fixture-provenance.json"
  )
  provenance <- jsonlite::fromJSON(provenance_path, simplifyVector = TRUE)
  if (!identical(as.integer(provenance$schema_version), 1L)) {
    stop("Unsupported cortical fixture provenance schema.", call. = FALSE)
  }
  resource <- provenance$resources[
    as.integer(provenance$resources$landmark_count) == landmark_count,
    , drop = FALSE
  ]
  if (nrow(resource) != 1L) {
    stop("Cortical fixture provenance has no unique requested size.",
         call. = FALSE)
  }
  landmark_path <- moment_validation_resource(resource$landmark_resource[[1L]])
  geodesic_path <- moment_validation_resource(resource$geodesic_resource[[1L]])
  observed_md5 <- unname(tools::md5sum(c(landmark_path, geodesic_path)))
  expected_md5 <- c(
    resource$landmark_resource_md5[[1L]],
    resource$geodesic_resource_md5[[1L]]
  )
  if (!identical(observed_md5, expected_md5)) {
    stop("Retained cortical fixture resource checksum mismatch.", call. = FALSE)
  }
  landmark_data <- utils::read.csv(
    landmark_path, stringsAsFactors = FALSE, check.names = FALSE
  )
  geodesic_data <- utils::read.csv(
    geodesic_path, stringsAsFactors = FALSE, check.names = FALSE
  )
  expected_index <- seq_len(landmark_count)
  if (!identical(as.integer(landmark_data$landmark_index), expected_index) ||
      !identical(as.integer(geodesic_data$landmark_index), expected_index) ||
      nrow(landmark_data) != landmark_count ||
      ncol(geodesic_data) != landmark_count + 1L) {
    stop("Retained cortical fixture dimensions or indices are invalid.",
         call. = FALSE)
  }
  geodesic <- as.matrix(geodesic_data[, -1L, drop = FALSE])
  storage.mode(geodesic) <- "double"
  vertices <- as.matrix(landmark_data[, c("x", "y", "z"), drop = FALSE])
  storage.mode(vertices) <- "double"
  weights <- as.numeric(landmark_data$weight)
  fixture <- list(
    vertices = vertices,
    weights = weights,
    geodesic = geodesic,
    landmark_vertex_ids = as.integer(landmark_data$full_vertex_id),
    dataset = "neuroatlas::fsaverage6",
    hemisphere = "left",
    surface = "pial",
    dataset_package_version = provenance$source$package_version,
    full_vertex_count = as.integer(provenance$source$full_vertex_count),
    full_face_count = as.integer(provenance$source$full_face_count),
    coordinate_scale_mm = provenance$source$coordinate_scale_mm,
    landmark_method = "farthest_point_on_inflated_surface",
    area_weight_method =
      "full_pial_face_area_aggregated_to_nearest_inflated_landmark",
    fixture_provenance_resource = basename(provenance_path),
    fixture_resource_md5 = observed_md5
  )
  fixture$fixture_fingerprint_md5 <- .moment_fixture_fingerprint(fixture)
  if (!identical(
      fixture$fixture_fingerprint_md5,
      resource$fixture_fingerprint_md5[[1L]]
    )) {
    stop("Retained cortical fixture semantic fingerprint mismatch.",
         call. = FALSE)
  }
  assign(cache_key, fixture, envir = .moment_validation_cache)
  fixture
}

.moment_surface_channels <- function(geodesic) {
  n <- nrow(geodesic)
  anchors <- unique(as.integer(round(seq(1, n, length.out = 8L))))
  if (length(anchors) < 7L) {
    stop("Cortical functional-map fixture needs at least seven anchors.",
         call. = FALSE)
  }
  distance_scale <- stats::median(geodesic[upper.tri(geodesic)])
  kernels <- vapply(anchors, function(anchor) {
    exp(-(geodesic[, anchor] / (.4 * distance_scale))^2)
  }, numeric(n))
  list(
    train = kernels[, c(1L, 3L, 5L, 7L), drop = FALSE],
    heldout = kernels[, c(2L, 4L, 6L), drop = FALSE],
    train_anchor_ids = anchors[c(1L, 3L, 5L, 7L)],
    heldout_anchor_ids = anchors[c(2L, 4L, 6L)]
  )
}

moment_make_surface_case <- function(seed, landmark_count = 36L) {
  set.seed(seed)
  mesh <- .moment_cortical_fixture(landmark_count)
  channels <- .moment_surface_channels(mesh$geodesic)
  n <- nrow(mesh$vertices)
  target_order <- sample(seq_len(n))
  angle <- runif(1L, -.7, .7)
  rotation <- matrix(
    c(cos(angle), -sin(angle), 0,
      sin(angle), cos(angle), 0,
      0, 0, 1),
    3L, 3L, byrow = TRUE
  )
  source_features <- channels$train + matrix(
    rnorm(length(channels$train), sd = .01), n, ncol(channels$train)
  )
  target_features <- channels$train[target_order, , drop = FALSE] + matrix(
    rnorm(length(channels$train), sd = .01), n, ncol(channels$train)
  )
  target_heldout <- channels$heldout[target_order, , drop = FALSE] + matrix(
    rnorm(length(channels$heldout), sd = .002), n, ncol(channels$heldout)
  )
  list(
    seed = seed,
    source_vertices = mesh$vertices,
    target_vertices = mesh$vertices[target_order, , drop = FALSE] %*% t(rotation),
    source_weights = mesh$weights,
    target_weights = mesh$weights[target_order],
    source_features = source_features,
    target_features = target_features,
    source_heldout = channels$heldout,
    target_heldout = target_heldout,
    source_geodesic = mesh$geodesic,
    target_geodesic = mesh$geodesic[target_order, target_order, drop = FALSE],
    truth_target_index = match(seq_len(n), target_order),
    target_order = target_order,
    area_weight_method = mesh$area_weight_method,
    dataset = mesh$dataset,
    hemisphere = mesh$hemisphere,
    surface = mesh$surface,
    dataset_package_version = mesh$dataset_package_version,
    fixture_fingerprint_md5 = mesh$fixture_fingerprint_md5,
    fixture_provenance_resource = mesh$fixture_provenance_resource,
    fixture_resource_md5 = mesh$fixture_resource_md5,
    full_vertex_count = mesh$full_vertex_count,
    full_face_count = mesh$full_face_count,
    coordinate_scale_mm = mesh$coordinate_scale_mm,
    landmark_method = mesh$landmark_method,
    landmark_vertex_ids = mesh$landmark_vertex_ids,
    train_anchor_ids = channels$train_anchor_ids,
    heldout_anchor_ids = channels$heldout_anchor_ids,
    functional_map_source = "frozen_geodesic_radial_basis_channels"
  )
}

.moment_surface_landmark_count <- function(protocol) {
  if (!is.null(protocol$landmark_count)) {
    return(as.integer(protocol$landmark_count))
  }
  as.integer(protocol$mesh_side)^2L
}

.moment_validate_surface_fixture <- function(problem, protocol) {
  checks <- c(
    dataset = identical(problem$dataset, protocol$dataset),
    package_version = identical(
      problem$dataset_package_version, protocol$dataset_package_version
    ),
    fixture_fingerprint = identical(
      problem$fixture_fingerprint_md5, protocol$fixture_fingerprint_md5
    ),
    fixture_provenance_resource = identical(
      problem$fixture_provenance_resource,
      protocol$fixture_provenance_resource
    ),
    hemisphere = identical(problem$hemisphere, protocol$hemisphere),
    surface = identical(problem$surface, protocol$surface),
    full_vertex_count = problem$full_vertex_count == protocol$full_vertex_count,
    full_face_count = problem$full_face_count == protocol$full_face_count,
    landmark_count = nrow(problem$source_vertices) == protocol$landmark_count,
    landmark_method = identical(
      problem$landmark_method, protocol$landmark_method
    ),
    area_weight_method = identical(
      problem$area_weight_method, protocol$area_weight_method
    ),
    functional_map_source = identical(
      problem$functional_map_source, protocol$functional_map_source
    ),
    train_anchor_ids = identical(
      as.integer(problem$train_anchor_ids), as.integer(protocol$train_anchor_ids)
    ),
    heldout_anchor_ids = identical(
      as.integer(problem$heldout_anchor_ids),
      as.integer(protocol$heldout_anchor_ids)
    )
  )
  if (!all(checks)) {
    stop(
      "Frozen cortical fixture mismatch: ",
      paste(names(checks)[!checks], collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.moment_surface_metrics <- function(fit, problem, geometry_kind, embedding_rank) {
  applied <- .moment_normalized_apply(fit$plans$sample, problem$target_heldout)
  active <- applied$active
  weights <- problem$source_weights[active]
  plan <- transport_plan_materialize(fit$plans$sample)
  row_mass <- rowSums(plan)
  geodesic_error <- numeric(nrow(plan))
  for (source in seq_len(nrow(plan))) {
    if (row_mass[[source]] > 0) {
      geodesic_error[[source]] <- sum(
        plan[source, ] * problem$target_geodesic[
          , problem$truth_target_index[[source]]
        ]
      ) / row_mass[[source]]
    } else {
      geodesic_error[[source]] <- Inf
    }
  }
  target_mass <- colSums(plan)
  discrepancy <- fit$coupling_discrepancy
  data.frame(
    solver_success = isTRUE(fit$final_stationarity_certified) &&
      isTRUE(fit$inner_uot_certified),
    status = fit$status,
    failure_reason = paste(fit$certification_failures, collapse = ";"),
    geometry_kind = geometry_kind,
    embedding_rank = embedding_rank,
    geometry_certified = isTRUE(fit$geometry_certified),
    geometry_exact = isTRUE(fit$geometry_exact),
    geometry_relative_error = fit$geometry_relative_error,
    geometry_heldout_error = if (isTRUE(fit$geometry_exact)) {
      0
    } else fit$geometry_heldout_error,
    heldout_correlation = .moment_weighted_correlation(
      applied$prediction[active, , drop = FALSE],
      problem$source_heldout[active, , drop = FALSE],
      weights
    ),
    heldout_r2 = .moment_weighted_r2(
      applied$prediction[active, , drop = FALSE],
      problem$source_heldout[active, , drop = FALSE],
      weights
    ),
    expected_geodesic_endpoint_error = .moment_weighted_mean(
      geodesic_error[active], weights
    ),
    transported_mass = sum(row_mass),
    discarded_source_mass = max(0, 1 - sum(row_mass)),
    discarded_target_mass = max(0, 1 - sum(target_mass)),
    coupling_l1_discrepancy = discrepancy$l1,
    coupling_max_discrepancy = discrepancy$max,
    validation_materialization = TRUE,
    dataset = problem$dataset,
    hemisphere = problem$hemisphere,
    surface = problem$surface,
    dataset_package_version = problem$dataset_package_version,
    fixture_fingerprint_md5 = problem$fixture_fingerprint_md5,
    full_vertex_count = problem$full_vertex_count,
    full_face_count = problem$full_face_count,
    coordinate_scale_mm = problem$coordinate_scale_mm,
    landmark_method = problem$landmark_method,
    landmark_vertex_ids = paste(problem$landmark_vertex_ids, collapse = ";"),
    area_weight_method = problem$area_weight_method,
    functional_map_source = problem$functional_map_source,
    train_anchor_ids = paste(problem$train_anchor_ids, collapse = ";"),
    heldout_anchor_ids = paste(problem$heldout_anchor_ids, collapse = ";"),
    stringsAsFactors = FALSE
  )
}

.moment_fit_surface <- function(problem, Cx, Cy, protocol) {
  fugw_factorized(
    Cx, Cy,
    wx = problem$source_weights,
    wy = problem$target_weights,
    M = sqeuclidean_cost(problem$source_features, problem$target_features),
    feature_weight = protocol$feature_weight,
    structure_weight = protocol$structure_weight,
    reg_marginals = protocol$reg_marginals,
    epsilon = protocol$epsilon,
    max_iter = as.integer(protocol$max_iter),
    tol = protocol$tol,
    max_iter_ot = as.integer(protocol$max_iter_ot),
    tol_ot = protocol$tol_ot,
    block_size = as.integer(protocol$block_size)
  )
}

moment_run_surface_case <- function(seed, protocol = NULL) {
  if (is.null(protocol)) protocol <- moment_validation_protocol()$surface
  problem <- moment_make_surface_case(
    seed, .moment_surface_landmark_count(protocol)
  )
  .moment_validate_surface_fixture(problem, protocol)
  tryCatch({
    fit <- .moment_fit_surface(
      problem,
      sqeuclidean_cost(problem$source_vertices),
      sqeuclidean_cost(problem$target_vertices),
      protocol
    )
    row <- .moment_surface_metrics(fit, problem, "exact_sqeuclidean", 3L)
    row$seed <- seed
    row$n_source <- nrow(problem$source_vertices)
    row$n_target <- nrow(problem$target_vertices)
    row
  }, error = function(error) {
    data.frame(
      solver_success = FALSE,
      status = "error",
      failure_reason = conditionMessage(error),
      geometry_kind = "exact_sqeuclidean",
      embedding_rank = 3L,
      geometry_certified = FALSE,
      geometry_exact = TRUE,
      geometry_relative_error = 0,
      geometry_heldout_error = 0,
      heldout_correlation = NA_real_,
      heldout_r2 = NA_real_,
      expected_geodesic_endpoint_error = NA_real_,
      transported_mass = NA_real_,
      discarded_source_mass = NA_real_,
      discarded_target_mass = NA_real_,
      coupling_l1_discrepancy = NA_real_,
      coupling_max_discrepancy = NA_real_,
      validation_materialization = FALSE,
      dataset = problem$dataset,
      hemisphere = problem$hemisphere,
      surface = problem$surface,
      dataset_package_version = problem$dataset_package_version,
      fixture_fingerprint_md5 = problem$fixture_fingerprint_md5,
      full_vertex_count = problem$full_vertex_count,
      full_face_count = problem$full_face_count,
      coordinate_scale_mm = problem$coordinate_scale_mm,
      landmark_method = problem$landmark_method,
      landmark_vertex_ids = paste(problem$landmark_vertex_ids, collapse = ";"),
      area_weight_method = problem$area_weight_method,
      functional_map_source = problem$functional_map_source,
      train_anchor_ids = paste(problem$train_anchor_ids, collapse = ";"),
      heldout_anchor_ids = paste(problem$heldout_anchor_ids, collapse = ";"),
      seed = seed,
      n_source = nrow(problem$source_vertices),
      n_target = nrow(problem$target_vertices),
      stringsAsFactors = FALSE
    )
  })
}

.moment_classical_embedding <- function(squared_distance, rank) {
  n <- nrow(squared_distance)
  center <- diag(n) - matrix(1 / n, n, n)
  gram <- -.5 * center %*% squared_distance %*% center
  decomposition <- eigen(gram, symmetric = TRUE)
  rank <- min(as.integer(rank), n - 1L)
  values <- pmax(decomposition$values[seq_len(rank)], 0)
  sweep(
    decomposition$vectors[, seq_len(rank), drop = FALSE],
    2L,
    sqrt(values),
    "*"
  )
}

moment_run_surface_rank_curve <- function(protocol = NULL) {
  if (is.null(protocol)) protocol <- moment_validation_protocol()$surface
  curve <- protocol$embedding_curve
  problem <- moment_make_surface_case(
    as.integer(curve$seed), .moment_surface_landmark_count(protocol)
  )
  .moment_validate_surface_fixture(problem, protocol)
  source_reference <- problem$source_geodesic^2
  target_reference <- problem$target_geodesic^2
  rows <- lapply(curve$ranks, function(rank) {
    source_embedding <- .moment_classical_embedding(source_reference, rank)
    target_embedding <- source_embedding[
      problem$target_order, , drop = FALSE
    ]
    Cx <- embedded_geometry_cost(
      source_embedding,
      source_reference,
      method = "classical_mds_mesh_geodesic",
      reference_metric = "squared_mesh_shortest_path_geodesic",
      seed = as.integer(curve$seed),
      relative_error_tolerance = curve$relative_error_tolerance,
      provenance = list(requested_rank = rank)
    )
    Cy <- embedded_geometry_cost(
      target_embedding,
      target_reference,
      method = "classical_mds_mesh_geodesic",
      reference_metric = "squared_mesh_shortest_path_geodesic",
      seed = as.integer(curve$seed) + 1L,
      relative_error_tolerance = curve$relative_error_tolerance,
      provenance = list(requested_rank = rank)
    )
    tryCatch({
      fit <- .moment_fit_surface(problem, Cx, Cy, protocol)
      row <- .moment_surface_metrics(
        fit, problem, "approximate_mesh_geodesic_mds", as.integer(rank)
      )
      row$seed <- as.integer(curve$seed)
      row$n_source <- nrow(problem$source_vertices)
      row$n_target <- nrow(problem$target_vertices)
      row
    }, error = function(error) {
      data.frame(
        solver_success = FALSE,
        status = "error",
        failure_reason = conditionMessage(error),
        geometry_kind = "approximate_mesh_geodesic_mds",
        embedding_rank = as.integer(rank),
        geometry_certified = FALSE,
        geometry_exact = FALSE,
        geometry_relative_error = NA_real_,
        geometry_heldout_error = NA_real_,
        heldout_correlation = NA_real_,
        heldout_r2 = NA_real_,
        expected_geodesic_endpoint_error = NA_real_,
        transported_mass = NA_real_,
        discarded_source_mass = NA_real_,
        discarded_target_mass = NA_real_,
        coupling_l1_discrepancy = NA_real_,
        coupling_max_discrepancy = NA_real_,
        validation_materialization = FALSE,
        dataset = problem$dataset,
        hemisphere = problem$hemisphere,
        surface = problem$surface,
        dataset_package_version = problem$dataset_package_version,
        fixture_fingerprint_md5 = problem$fixture_fingerprint_md5,
        full_vertex_count = problem$full_vertex_count,
        full_face_count = problem$full_face_count,
        coordinate_scale_mm = problem$coordinate_scale_mm,
        landmark_method = problem$landmark_method,
        landmark_vertex_ids = paste(
          problem$landmark_vertex_ids, collapse = ";"
        ),
        area_weight_method = problem$area_weight_method,
        functional_map_source = problem$functional_map_source,
        train_anchor_ids = paste(problem$train_anchor_ids, collapse = ";"),
        heldout_anchor_ids = paste(problem$heldout_anchor_ids, collapse = ";"),
        seed = as.integer(curve$seed),
        n_source = nrow(problem$source_vertices),
        n_target = nrow(problem$target_vertices),
        stringsAsFactors = FALSE
      )
    })
  })
  do.call(rbind, rows)
}

moment_surface_dense_parity <- function(protocol = NULL) {
  if (is.null(protocol)) protocol <- moment_validation_protocol()$surface
  problem <- moment_make_surface_case(
    7101L, landmark_count = as.integer(protocol$dense_parity_landmark_count)
  )
  fixture_fingerprint_passed <- identical(
    problem$fixture_fingerprint_md5,
    protocol$dense_parity_fixture_fingerprint_md5
  )
  Cx <- sqeuclidean_cost(problem$source_vertices)
  Cy <- sqeuclidean_cost(problem$target_vertices)
  M <- sqeuclidean_cost(problem$source_features, problem$target_features)
  factorized <- .moment_fit_surface(problem, Cx, Cy, protocol)
  dense <- fugw_kl(
    cost_block(Cx),
    cost_block(Cy),
    wx = problem$source_weights,
    wy = problem$target_weights,
    M = cost_block(M),
    feature_weight = protocol$feature_weight,
    structure_weight = protocol$structure_weight,
    reg_marginals = protocol$reg_marginals,
    epsilon = protocol$epsilon,
    max_iter = as.integer(protocol$max_iter),
    tol = protocol$tol,
    max_iter_ot = as.integer(protocol$max_iter_ot),
    tol_ot = protocol$tol_ot,
    precision = "strict_double"
  )
  probe <- cbind(
    sin(seq_len(nrow(problem$target_vertices))),
    cos(seq_len(nrow(problem$target_vertices))),
    seq_len(nrow(problem$target_vertices)) /
      nrow(problem$target_vertices)
  )
  action_error <- max(abs(
    transport_plan_apply(factorized$plans$sample, probe) -
      dense$pi_samp %*% probe
  ))
  objective_error <- abs(factorized$fugw_cost - dense$fugw_cost)
  structure_error <- abs(
    factorized$objective_decomposition$structure_unweighted -
      dense$objective_decomposition$structure_unweighted
  )
  feature_error <- abs(
    factorized$objective_decomposition$feature_unweighted -
      dense$objective_decomposition$feature_unweighted
  )
  regularization_error <- abs(
    factorized$objective_decomposition$regularization -
      dense$objective_decomposition$regularization
  )
  maximum_error <- max(
    objective_error, structure_error, feature_error,
    regularization_error, action_error
  )
  data.frame(
    seed = 7101L,
    n_source = nrow(problem$source_vertices),
    fixture_fingerprint_md5 = problem$fixture_fingerprint_md5,
    fixture_fingerprint_passed = fixture_fingerprint_passed,
    factorized_certified = isTRUE(factorized$final_stationarity_certified),
    dense_converged = isTRUE(dense$converged),
    objective_absolute_error = objective_error,
    structure_component_error = structure_error,
    feature_component_error = feature_error,
    regularization_component_error = regularization_error,
    sample_plan_action_max_error = action_error,
    maximum_parity_error = maximum_error,
    tolerance = protocol$dense_parity_tolerance,
    passed = fixture_fingerprint_passed &&
      isTRUE(factorized$final_stationarity_certified) &&
      isTRUE(dense$converged) &&
      maximum_error <= protocol$dense_parity_tolerance,
    stringsAsFactors = FALSE
  )
}

moment_run_surface_suite <- function(
    seeds = NULL,
    output = file.path("inst", "bench", "moment-fugw-surface-evidence.csv"),
    rank_output = file.path(
      "inst", "bench", "moment-fugw-geometry-rank-evidence.csv"
    )) {
  protocol <- moment_validation_protocol()$surface
  if (is.null(seeds)) seeds <- protocol$evaluation_seeds
  rows <- do.call(rbind, lapply(seeds, function(seed) {
    message("Surface validation seed ", seed)
    moment_run_surface_case(seed, protocol)
  }))
  rows$evaluation_split <- ifelse(
    rows$seed %in% protocol$evaluation_seeds, "frozen_evaluation", "tuning"
  )
  parity <- moment_surface_dense_parity(protocol)
  rows$dense_parity_passed <- isTRUE(parity$passed)
  rows$dense_parity_maximum_error <- parity$maximum_parity_error
  rows$dense_parity_tolerance <- parity$tolerance
  rows$suite_gate <- all(rows$solver_success) && isTRUE(parity$passed) &&
    all(rows$geometry_certified) && all(rows$geometry_exact) &&
    all(is.finite(rows$heldout_correlation)) &&
    all(is.finite(rows$expected_geodesic_endpoint_error))
  rank_rows <- moment_run_surface_rank_curve(protocol)
  rank_rows$numerical_and_geometry_error_separate <- TRUE
  full_protocol <- moment_validation_protocol()
  rows$protocol_schema_version <- full_protocol$schema_version
  rows$protocol_status <- full_protocol$status
  rank_rows$protocol_schema_version <- full_protocol$schema_version
  rank_rows$protocol_status <- full_protocol$status
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(rows, output, row.names = FALSE)
  utils::write.csv(rank_rows, rank_output, row.names = FALSE)
  attr(rows, "dense_parity") <- parity
  attr(rows, "rank_curve") <- rank_rows
  rows
}
