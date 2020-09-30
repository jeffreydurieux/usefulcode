# Thu Aug  8 10:28:46 2019
# Author: Jeffrey Durieux, MSc
# Still work in progress
# What: make a function that customizes the linting process of my code
library(lintr)

source_file <- get_source_expressions('~/Repositories/CodeFromTom/CLPAR_CheckEmptyClusters (1).R')

### available linters ###  
absolute_path_linter(source_file)

assignment_linter(source_file)

closed_curly_linter(allow_single_line = FALSE)

commas_linter(source_file)

commented_code_linter(source_file)

infix_spaces_linter(source_file)

line_length_linter(length)

no_tab_linter(source_file)

object_usage_linter(source_file)

camel_case_linter(source_file)

snake_case_linter(source_file)
