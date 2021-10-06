#!/usr/bin/Rscript

# Function: plus two variable

## Usage: Rscript R_script_template.R \
#         --first_number 1 \ 
#         --second_number 2

# optional arguments:
#    -h, --help  show this message and exit
#    -f, --first_number Input first number
#    -s, --second_number Input second number

# R package
library("argparser")

# Create a parser
p <- arg_parser("plus two number")

# Add command line arguments
p <- add_argument(p, "--first_number",  help="first number to plus", type="numeric", default = 2)
p <- add_argument(p, "--second_number", help="second number to plus", type="numeric", default = 3)

# Parse the command line arguments
argv <- parse_args(p)

# Do work based on the passed arguments
cat(paste0("First number: ", argv$first_number, "\n"))
cat(paste0("Second number: ", argv$second_number, "\n"))
cat(paste0("Total number: ",argv$first_number + argv$second_number, "\n"))