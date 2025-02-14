
# conda is causing wierd problems; check always abort R: https://forum.posit.co/t/unable-to-install-build-tools-when-doing-devtools-check/183308
# run `conda config --set auto_activate_base false` in command to prevent base conda to be activated.
# second debug working well -- `mv ~/.R/Makevars ~/.R/backup_Makevars`; it seems R is pointing to the wrong clang location using the Makevars file; this
# is connectd to the conda env issue, which is the reason fixing it first is helpful; but later they learned to do it wrong again by the Makevars.
# worth checking `which clang` in the command line; here it will say `/Users/xilin/anaconda3/bin/clang` which is not good...

# load dependency
usethis::use_package_doc() # this needs to start first (otherwise would remove all documents)
usethis::use_package("dplyr", min_version = "1.1.0")
usethis::use_import_from("dplyr", c("anti_join", "arrange", "group_by",
                                    "bind_rows", "filter", "left_join",
                                    "mutate", "n", "pull", "rename",
                                    "select", "slice", "summarise", "ungroup",
                                    "data_frame", "if_else", "add_row"))
usethis::use_import_from("stats", c("cov", "sd", "var"))
# usethis::use_package("data.table")
# usethis::use_import_from("data.table", c("data.table"))
usethis::use_pipe(export = T)
usethis::use_tidy_description()
usethis::use_gpl_license() # force other package use my code to be open-source as well.

# create data set creation documents -- should move to data-raw folder to edit how data is generated
usethis::use_data_raw("create_nonidentify_data")


# create tests
usethis::use_testthat(3)
usethis::use_test("BPR_model_evaluation")
