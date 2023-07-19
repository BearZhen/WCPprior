test_that("testing arc length function", {
  theta <- seq(from=pi,to=3*pi/2,length.out = 50)

  coords <- cbind(sin(theta),1+ cos(theta))
  
  arclengths <- compute_partial_arc_lengths(coords)[,3]

  true_arclengths <- (theta-pi)
  error_arc <- sum((true_arclengths - arclengths)^2)
  testthat::expect_equal(error_arc, 0, tolerance = 1e-7)
})
