rfugw_valid_commit <- function(value) {
  is.character(value) && length(value) == 1L && !is.na(value) &&
    grepl("^[[:xdigit:]]{40}$", value)
}

rfugw_normalize_commit <- function(value, label = "source commit") {
  if (!rfugw_valid_commit(value)) {
    stop(label, " must be exactly 40 hexadecimal characters.", call. = FALSE)
  }
  tolower(value)
}

rfugw_git_provenance <- function(repository_root = ".") {
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
      rfugw_valid_commit(commit_result$value[[1L]])) {
    tolower(commit_result$value[[1L]])
  } else {
    NA_character_
  }
  status_ok <- isTRUE(status_result$ok)
  complete <- !is.na(commit) && status_ok
  dirty <- if (status_ok) length(status_result$value) > 0L else NA
  list(
    repository_root = repository_root,
    commit = commit,
    working_tree_dirty = dirty,
    git_status_entry_count = if (status_ok) {
      as.integer(length(status_result$value))
    } else {
      NA_integer_
    },
    git_provenance_complete = complete,
    exact_commit_evidence = complete && identical(dirty, FALSE),
    git_commit_command_status = commit_result$status,
    git_status_command_status = status_result$status
  )
}
