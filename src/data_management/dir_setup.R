##Directory setup:
#Aim: Create a subdirectory structure to run snakemake without the need to
# change the config file

args <- (commandArgs(trailingOnly = TRUE))
for (i in seq_len(length(args))) {
  eval(parse(text = args[[i]]))
}

## Show arguments
print(dir_names)
dir_names <- strsplit(gsub(" ","",dir_names), ",")[[1]]
print(dir_names)
print(class(dir_names))


if( !dir.exists("out") ){
  stop("You can only automtically build subdirectories, if a `log` and an
       `output` directory already exist")
}

if( !dir.exists("log") ){
  stop("You can only automtically build subdirectories, if a `log` and an
       `output` directory already exist")
}

lapply(dir_names, function(name){
  #create output directories
  dir.create(paste0("out/", name))
  #create log directories
  dir.create(paste0("log/", name))
})
