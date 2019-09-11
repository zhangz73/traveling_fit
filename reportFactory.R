library(reportfactory)

destination <- file.path(getwd(), "my_reports")
new_factory(destination)

update_reports()
list_outputs()

list_reports()

compile_report("models_areaId2ad4_2019-08-28.Rmd", quiet = T)

compile_report("pfsi_sim_tools_2019-09-11.Rmd", quiet = T)
