# +------------------------------------------------------+
# | Utility Functions                                    |
# +------------------------------------------------------+

fls = list.files("./C_cea_code/code00c_functions/", pattern="^[fcn]")
for (i in 1:length(fls)){source(paste0("./C_cea_code/code00c_functions/", fls[i]))}
