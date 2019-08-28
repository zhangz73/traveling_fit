library(reportfactory)

destination <- file.path(getwd(), "my_reports")
new_factory(destination)

update_reports()
list_outputs()

list_reports()

compile_report("models_areaId2ad4_2019-08-28.Rmd", quiet = T)
