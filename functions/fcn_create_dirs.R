create_dir_if_necessary = function(x){
  # This function checks whether a path directory exists (defined by string x) and creates it if not.
  # It supports directory trees of multiple levels: if any folder does not exist, 
  # then it starts building the directory tree to create the folder described by x. 
  # However, all folders that exist (and their contents) are undisturbed by this procedure.
  
  tmp1 = strsplit(x, "/")[[1]]
  tmp1 = tmp1[tmp1!=""]
  for (i in 1:length(tmp1)){
    tmp2 = paste(tmp1[1:i], sep="", collapse="/")
    # print(tmp2)
    if(!dir.exists(file.path(tmp2))){dir.create(file.path(tmp2))}
  }

  # old, simple version: function(x){if(!dir.exists(file.path(x))){dir.create(file.path(x))}}
  
}
