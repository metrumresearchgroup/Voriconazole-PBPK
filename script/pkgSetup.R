## This script installs all required packages
## Before you proceed make the script directory your working directory ith the function setwd()
scriptDir <- getwd()
pkgDir <- file.path(scriptDir, "pkg")
libDir <- file.path(scriptDir, "lib")
dir.create(libDir)

pkgs <- list.files(pkgDir)
pkgss <- c(paste(pkgDir, "/", pkgs, sep=""))

install.packages(pkgss, repos=NULL, type="source", lib=libDir)

