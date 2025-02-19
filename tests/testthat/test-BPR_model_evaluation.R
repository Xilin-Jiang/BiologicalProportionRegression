test_that("wrapper function with multiple ", {
  target_traits <- cbind(X_EXAMPLE,X_EXAMPLE, X_EXAMPLE)
  helper_traits <- cbind(Y_EXAMPLE,Y_EXAMPLE, Z_EXAMPLE[,1])
  auxiliary_traits <- Z_EXAMPLE
  BPR_fast_wrapper(target_traits, helper_traits, auxiliary_traits, bootstrapping_number = 5)
})
