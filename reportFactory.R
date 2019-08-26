library(reportfactory)

destination <- file.path(getwd(), "my_reports")
new_factory(destination)

update_reports()
list_outputs()

compile_report("models_areaId2ad4_2019-08-26", quiet = T)
compile_report("useful_functions_summary-2019-08-10", quiet = T)