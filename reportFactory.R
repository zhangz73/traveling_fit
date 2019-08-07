library(reportfactory)

destination <- file.path(getwd(), "my_reports")
new_factory(destination)

update_reports()
list_outputs()