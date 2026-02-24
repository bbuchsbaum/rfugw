test_that("non-entropic semirelaxed matches POT fixture", {
  fx <- read_fixture("pot_api_ports_fixture.json")
  fxs <- fx$semirelaxed

  out_srgw <- rfugw::semirelaxed_gromov_wasserstein(
    C1 = fxs$inputs$C1,
    C2 = fxs$inputs$C2,
    p = fxs$inputs$p,
    max_iter = fxs$params$max_iter,
    tol_rel = fxs$params$tol_rel,
    tol_abs = fxs$params$tol_abs
  )
  out_srfgw <- rfugw::semirelaxed_fused_gromov_wasserstein(
    M = fxs$inputs$M,
    C1 = fxs$inputs$C1,
    C2 = fxs$inputs$C2,
    p = fxs$inputs$p,
    alpha = fxs$params$alpha,
    max_iter = fxs$params$max_iter,
    tol_rel = fxs$params$tol_rel,
    tol_abs = fxs$params$tol_abs
  )

  expect_equal(out_srgw$srgw_dist, fxs$outputs$srgw_dist, tolerance = 5e-4)
  expect_equal(out_srfgw$srfgw_dist, fxs$outputs$srfgw_dist, tolerance = 5e-4)
  expect_equal(unname(out_srgw$plan), unname(fxs$outputs$T_srgw), tolerance = 3e-3)
  expect_equal(unname(out_srfgw$plan), unname(fxs$outputs$T_srfgw), tolerance = 3e-3)
})

test_that("entropic partial GW/FGW matches POT fixture", {
  fx <- read_fixture("pot_api_ports_fixture.json")
  fxe <- fx$entropic_partial

  out_gw <- rfugw::entropic_partial_gromov_wasserstein(
    C1 = fxe$inputs$C1,
    C2 = fxe$inputs$C2,
    p = fxe$inputs$p,
    q = fxe$inputs$q,
    reg = fxe$params$reg,
    m = fxe$params$m,
    numItermax = fxe$params$numItermax,
    tol = fxe$params$tol,
    log = TRUE
  )
  out_fgw <- rfugw::entropic_partial_fused_gromov_wasserstein(
    M = fxe$inputs$M,
    C1 = fxe$inputs$C1,
    C2 = fxe$inputs$C2,
    p = fxe$inputs$p,
    q = fxe$inputs$q,
    reg = fxe$params$reg,
    m = fxe$params$m,
    alpha = fxe$params$alpha,
    numItermax = fxe$params$numItermax,
    tol = fxe$params$tol,
    log = TRUE
  )

  expect_lt(abs(out_gw$partial_gw_dist - fxe$outputs$epgw_dist), 3e-2)
  expect_lt(abs(out_fgw$partial_fgw_dist - fxe$outputs$epfgw_dist), 3e-2)
  expect_equal(unname(out_gw$plan), unname(fxe$outputs$T_epgw), tolerance = 2e-2)
  expect_equal(unname(out_fgw$plan), unname(fxe$outputs$T_epfgw), tolerance = 2e-2)
})

test_that("UCOOT KL sinkhorn matches POT fixture", {
  fx <- read_fixture("pot_api_ports_fixture.json")
  fxu <- fx$ucoot

  out <- rfugw::unbalanced_co_optimal_transport(
    X = fxu$inputs$X,
    Y = fxu$inputs$Y,
    reg_marginals = unlist(fxu$params$reg_marginals),
    epsilon = unlist(fxu$params$epsilon),
    divergence = "kl",
    unbalanced_solver = "sinkhorn",
    max_iter = fxu$params$max_iter,
    tol = fxu$params$tol,
    max_iter_ot = fxu$params$max_iter_ot,
    tol_ot = fxu$params$tol_ot,
    log = TRUE
  )

  expect_equal(out$ucoot_cost, fxu$outputs$ucoot_cost, tolerance = 2e-2)
  expect_equal(unname(out$pi_samp), unname(fxu$outputs$pi_samp), tolerance = 2e-2)
  expect_equal(unname(out$pi_feat), unname(fxu$outputs$pi_feat), tolerance = 2e-2)
})

test_that("sampled GW fixture comparison is in-range", {
  fx <- read_fixture("pot_api_ports_fixture.json")
  fxs <- fx$sampled

  out <- rfugw::sampled_gromov_wasserstein(
    C1 = fxs$inputs$C1,
    C2 = fxs$inputs$C2,
    p = fxs$inputs$p,
    q = fxs$inputs$q,
    nb_samples_grad = unlist(fxs$params$nb_samples_grad),
    epsilon = fxs$params$epsilon,
    max_iter = fxs$params$max_iter,
    random_state = fxs$params$random_state,
    log = TRUE
  )

  expect_equal(out$gw_dist_estimated, fxs$outputs$gw_dist_estimated, tolerance = 2e-1)
  expect_equal(unname(out$plan), unname(fxs$outputs$plan), tolerance = 2e-1)
})
