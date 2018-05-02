## This script installs all required packages and creates a "lib" directory
## Before you proceed make the script directory your working directory with the function setwd()
scriptDir <- getwd()
pkgDir <- file.path(scriptDir, "pkg")
libDir <- file.path(scriptDir, "lib")
dir.create(libDir)

pkgs <- list.files(pkgDir)
pkgss <- c(paste(pkgDir, "/", pkgs, sep=""))

install.packages(pkgss, repos=NULL, type="source", lib=libDir)

