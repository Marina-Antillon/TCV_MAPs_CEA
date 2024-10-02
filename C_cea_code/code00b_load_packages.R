# +------------------------------------------------------+
# | Load packages                                        |
# +------------------------------------------------------+

# Some packages like to announce they're loaded, so
# suppressPackageStartupMessages gets them to be quiet while still letting
# warnings/errors through

# For full features of pipes
library(magrittr)

# For plotting
suppressPackageStartupMessages(library(ggplot2))

# For data manipulation
library(dplyr, warn.conflicts=FALSE) # Masks stats:filter,lag, base:intersect,setdiff,setequal,union
library(cubelyr)
library(tidyr)
library(abind)

# For loading data
library(readxl)
