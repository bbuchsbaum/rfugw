args <- commandArgs(trailingOnly = TRUE)
surface_dir <- if (length(args) >= 1L) {
  normalizePath(args[[1L]], mustWork = TRUE)
} else {
  stop(
    paste0(
      "Usage: Rscript inst/bench/generate_moment_fugw_cortical_fixture.R ",
      "/path/to/neuroatlas/data-raw/fsaverage6 [output_dir]"
    ),
    call. = FALSE
  )
}
output_dir <- if (length(args) >= 2L) {
  normalizePath(args[[2L]], mustWork = TRUE)
} else {
  normalizePath(file.path("inst", "bench"), mustWork = TRUE)
}

if (!requireNamespace("igraph", quietly = TRUE) ||
    !requireNamespace("jsonlite", quietly = TRUE)) {
  stop("Fixture generation requires igraph and jsonlite.", call. = FALSE)
}

validation <- new.env(parent = globalenv())
sys.source(
  file.path(output_dir, "moment_fugw_scientific_validation.R"),
  envir = validation
)

read_freesurfer_triangle_surface <- function(path) {
  connection <- file(path, open = "rb")
  on.exit(close(connection), add = TRUE)
  magic_raw <- readBin(connection, what = "raw", n = 3L)
  if (length(magic_raw) != 3L) {
    stop("Truncated FreeSurfer surface: ", path, call. = FALSE)
  }
  magic <- sum(as.integer(magic_raw) * c(65536L, 256L, 1L))
  if (magic != 16777214L) {
    stop("Only FreeSurfer triangle surfaces are supported: ", path,
         call. = FALSE)
  }
  read_binary_line <- function() {
    bytes <- raw()
    repeat {
      byte <- readBin(connection, what = "raw", n = 1L)
      if (!length(byte)) stop("Truncated FreeSurfer header.", call. = FALSE)
      if (as.integer(byte) == 10L) break
      bytes <- c(bytes, byte)
    }
    rawToChar(bytes)
  }
  comment <- read_binary_line()
  stamp <- read_binary_line()
  counts <- readBin(
    connection, what = integer(), n = 2L, size = 4L,
    signed = TRUE, endian = "big"
  )
  if (length(counts) != 2L || any(counts <= 0L)) {
    stop("Invalid FreeSurfer vertex/face counts: ", path, call. = FALSE)
  }
  vertices <- matrix(readBin(
    connection, what = numeric(), n = 3L * counts[[1L]], size = 4L,
    endian = "big"
  ), ncol = 3L, byrow = TRUE)
  faces <- matrix(readBin(
    connection, what = integer(), n = 3L * counts[[2L]], size = 4L,
    signed = TRUE, endian = "big"
  ), ncol = 3L, byrow = TRUE) + 1L
  if (nrow(vertices) != counts[[1L]] || nrow(faces) != counts[[2L]] ||
      any(faces < 1L) || any(faces > nrow(vertices))) {
    stop("Truncated or invalid FreeSurfer surface payload: ", path,
         call. = FALSE)
  }
  list(vertices = vertices, faces = faces, comment = comment, stamp = stamp)
}

write_numeric_csv <- function(data, path) {
  encoded <- lapply(data, function(column) {
    if (is.integer(column)) return(as.character(column))
    if (is.numeric(column)) {
      return(sprintf("%a", column))
    }
    as.character(column)
  })
  encoded <- as.data.frame(encoded, check.names = FALSE,
                           stringsAsFactors = FALSE)
  utils::write.table(
    encoded, path, sep = ",", row.names = FALSE, col.names = TRUE,
    quote = FALSE, na = ""
  )
}

pial_path <- file.path(surface_dir, "lh.pial")
inflated_path <- file.path(surface_dir, "lh.inflated")
if (!file.exists(pial_path) || !file.exists(inflated_path)) {
  stop("surface_dir must contain lh.pial and lh.inflated.", call. = FALSE)
}
pial <- read_freesurfer_triangle_surface(pial_path)
inflated <- read_freesurfer_triangle_surface(inflated_path)
if (!identical(pial$faces, inflated$faces)) {
  stop("The pial and inflated surface topology differs.", call. = FALSE)
}

triangle_areas <- validation$.moment_triangle_areas(
  pial$vertices, pial$faces
)
vertex_areas <- validation$.moment_weighted_tabulate(
  as.vector(pial$faces),
  rep(triangle_areas / 3, times = 3L),
  nbins = nrow(pial$vertices)
)
center <- colMeans(pial$vertices)
centered <- sweep(pial$vertices, 2L, center, "-")
coordinate_scale_mm <- sqrt(mean(rowSums(centered^2)))
pial_vertices <- centered / coordinate_scale_mm
inflated_vertices <- scale(inflated$vertices, center = TRUE, scale = FALSE)
inflated_vertices <- inflated_vertices /
  sqrt(mean(rowSums(inflated_vertices^2)))

all_edges <- rbind(
  pial$faces[, c(1L, 2L), drop = FALSE],
  pial$faces[, c(2L, 3L), drop = FALSE],
  pial$faces[, c(1L, 3L), drop = FALSE]
)
edge_low <- pmin(all_edges[, 1L], all_edges[, 2L])
edge_high <- pmax(all_edges[, 1L], all_edges[, 2L])
edge_key <- edge_low + (edge_high - 1) * nrow(pial_vertices)
keep <- !duplicated(edge_key)
edges <- cbind(edge_low[keep], edge_high[keep])
edge_weights <- sqrt(rowSums(
  (pial_vertices[edges[, 1L], , drop = FALSE] -
     pial_vertices[edges[, 2L], , drop = FALSE])^2
))
graph <- igraph::graph_from_edgelist(edges, directed = FALSE)
graph <- igraph::set_edge_attr(graph, "weight", value = edge_weights)

fixture_metadata <- list(
  dataset = "neuroatlas::fsaverage6",
  hemisphere = "left",
  surface = "pial",
  dataset_package_version = "0.1.0.9000",
  full_vertex_count = nrow(pial_vertices),
  full_face_count = nrow(pial$faces),
  coordinate_scale_mm = coordinate_scale_mm,
  landmark_method = "farthest_point_on_inflated_surface",
  area_weight_method =
    "full_pial_face_area_aggregated_to_nearest_inflated_landmark"
)

make_fixture <- function(landmark_count) {
  landmarks <- validation$.moment_farthest_landmarks(
    inflated_vertices, landmark_count
  )
  weights <- validation$.moment_aggregate_area_weights(
    inflated_vertices, landmarks, vertex_areas
  )
  geodesic <- igraph::distances(
    graph, v = landmarks, to = landmarks,
    weights = igraph::E(graph)$weight, algorithm = "dijkstra"
  )
  fixture <- c(list(
    vertices = pial_vertices[landmarks, , drop = FALSE],
    weights = weights,
    geodesic = geodesic,
    landmark_vertex_ids = landmarks
  ), fixture_metadata)
  fixture$fixture_fingerprint_md5 <-
    validation$.moment_fixture_fingerprint(fixture)
  fixture
}

write_fixture <- function(fixture) {
  count <- nrow(fixture$vertices)
  landmark_name <- paste0(
    "moment-fugw-fsaverage6-", count, "-landmarks.csv"
  )
  geodesic_name <- paste0(
    "moment-fugw-fsaverage6-", count, "-geodesic.csv"
  )
  landmarks <- data.frame(
    landmark_index = seq_len(count),
    full_vertex_id = fixture$landmark_vertex_ids,
    x = fixture$vertices[, 1L],
    y = fixture$vertices[, 2L],
    z = fixture$vertices[, 3L],
    weight = fixture$weights,
    stringsAsFactors = FALSE
  )
  geodesic <- as.data.frame(fixture$geodesic, check.names = FALSE)
  names(geodesic) <- paste0("d", seq_len(count))
  geodesic <- cbind(landmark_index = seq_len(count), geodesic)
  write_numeric_csv(landmarks, file.path(output_dir, landmark_name))
  write_numeric_csv(geodesic, file.path(output_dir, geodesic_name))
  retained_landmarks <- utils::read.csv(
    file.path(output_dir, landmark_name), check.names = FALSE,
    stringsAsFactors = FALSE
  )
  retained_geodesic <- utils::read.csv(
    file.path(output_dir, geodesic_name), check.names = FALSE,
    stringsAsFactors = FALSE
  )
  retained_vertices <- as.matrix(
    retained_landmarks[, c("x", "y", "z"), drop = FALSE]
  )
  storage.mode(retained_vertices) <- "double"
  retained_geodesic_matrix <- as.matrix(
    retained_geodesic[, -1L, drop = FALSE]
  )
  storage.mode(retained_geodesic_matrix) <- "double"
  retained_fixture <- c(list(
    vertices = retained_vertices,
    weights = as.numeric(retained_landmarks$weight),
    geodesic = retained_geodesic_matrix,
    landmark_vertex_ids = as.integer(retained_landmarks$full_vertex_id)
  ), fixture_metadata)
  retained_fingerprint <- validation$.moment_fixture_fingerprint(
    retained_fixture
  )
  if (!identical(retained_fingerprint, fixture$fixture_fingerprint_md5)) {
    stop(
      "CSV round-trip changed the semantic fixture fingerprint: expected ",
      fixture$fixture_fingerprint_md5, ", observed ", retained_fingerprint,
      "; vertices identical=", identical(
        retained_fixture$vertices, fixture$vertices
      ), ", max diff=", max(abs(
        retained_fixture$vertices - fixture$vertices
      )), "; weights identical=", identical(
        retained_fixture$weights, fixture$weights
      ), ", max diff=", max(abs(
        retained_fixture$weights - fixture$weights
      )), "; geodesic identical=", identical(
        retained_fixture$geodesic, fixture$geodesic
      ), ", max diff=", max(abs(
        retained_fixture$geodesic - fixture$geodesic
      )),
      call. = FALSE
    )
  }
  list(
    landmark_count = count,
    fixture_fingerprint_md5 = fixture$fixture_fingerprint_md5,
    landmark_resource = landmark_name,
    landmark_resource_md5 = unname(tools::md5sum(
      file.path(output_dir, landmark_name)
    )),
    geodesic_resource = geodesic_name,
    geodesic_resource_md5 = unname(tools::md5sum(
      file.path(output_dir, geodesic_name)
    ))
  )
}

fixtures <- lapply(c(16L, 36L), make_fixture)
resources <- lapply(fixtures, write_fixture)
source_root <- normalizePath(file.path(surface_dir, "..", ".."),
                             mustWork = TRUE)
git_commit <- suppressWarnings(system2(
  "git", c("-C", source_root, "rev-parse", "HEAD"),
  stdout = TRUE, stderr = FALSE
))
if (!length(git_commit) || !is.null(attr(git_commit, "status"))) {
  git_commit <- NA_character_
}
provenance <- list(
  schema_version = 1L,
  generated_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  derivation = paste0(
    "Deterministic numeric landmark derivatives of the left fsaverage6 pial ",
    "and inflated surfaces. Pial coordinates are centered and RMS-normalized; ",
    "landmarks are farthest-point samples on the normalized inflated surface; ",
    "weights aggregate full-mesh pial triangle area; geodesics are shortest ",
    "paths over normalized-pial mesh edges."
  ),
  source = list(
    package = "neuroatlas",
    package_version = "0.1.0.9000",
    repository_commit = if (length(git_commit)) git_commit[[1L]] else NA_character_,
    raw_pial_md5 = unname(tools::md5sum(pial_path)),
    raw_inflated_md5 = unname(tools::md5sum(inflated_path)),
    full_vertex_count = nrow(pial_vertices),
    full_face_count = nrow(pial$faces),
    coordinate_scale_mm = coordinate_scale_mm,
    upstream_terms_note = paste0(
      "Derived from FreeSurfer fsaverage6 geometry bundled by neuroatlas; ",
      "the neuroatlas package is MIT-licensed, while upstream surface terms ",
      "remain authoritative."
    )
  ),
  resources = resources
)
jsonlite::write_json(
  provenance,
  file.path(output_dir, "moment-fugw-fsaverage6-fixture-provenance.json"),
  auto_unbox = TRUE, pretty = TRUE, na = "null", digits = 17
)

cat("Generated Moment-FUGW cortical fixtures in", output_dir, "\n")
for (resource in resources) {
  cat(
    resource$landmark_count, "landmarks:",
    resource$fixture_fingerprint_md5, "\n"
  )
}
