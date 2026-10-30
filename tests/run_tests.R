# Run the regression tests from the repo root: Rscript tests/run_tests.R
for (f in list.files("app/R", pattern = "\\.R$", full.names = TRUE)) source(f)
library(testthat)
test_dir("tests/testthat", reporter = "summary", stop_on_failure = TRUE)
