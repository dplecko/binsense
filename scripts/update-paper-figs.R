
# Get project root
root <- rprojroot::find_root(rprojroot::has_dir(".git"))

# List of R scripts or Rmd files to run
script_files <- c("nonpar.Rmd", "par.Rmd", "fos-gs-cond.R", "init-sense.R")
sel <- c(1)
script_files <- script_files[sel]

# Run the scripts or Rmd files
for (file in script_files) {
  
  if(grepl("\\.R$", file)) {
    source(file.path(root, "scripts", file))
  } else if(grepl("\\.Rmd$", file)) {
    rmarkdown::render(file.path(root, "vignettes", file))
  }
}

# Commit and push changes to Git
system(paste(
  "cd ~/ICU/op/paper &&
  git add figures/ &&",
  paste('git commit -m "Update images -', Sys.Date(), '" &&'),
  "git pull && git push"
))
