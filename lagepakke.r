# oppstart, tilgang til pakker anbefalt av Gandrud
# og Wickham, og Rmpfr:

# bare kjøre install.packages en gang:
install.packages("brew","countrycode","devtools","dplyr","ggplot2",
                 "googleVis","knitr","MCMCpack","repmis","RCurl",
                 "rmarkdown","texreg","tidyr","WDI",
                 "xtable","Zelig", "roxygen2","Rmpfr")
# virker ikke, tuller med at installasjonsomr?det for R er
# skrivebeskyttet selv om jeg har valgt det i egen mappe og
# fjernet skrivebeskyttelsen, det blir dill, den g?r ikke
# vekk likevel. Ikke ok, jeg kj?rer dette som adm-bruker.

# prøver å følge et råd på
# https://support.rstudio.com/hc/en-us/community/posts/205370887-install-packages-problem

install.packages("brew","countrycode","devtools","dplyr","ggplot2",
                 "googleVis","knitr","MCMCpack","repmis","RCurl",
                 "rmarkdown","texreg","tidyr","WDI",
                 "xtable","Zelig", "roxygen2","Rmpfr",
                 repos = getOption("repos")[["CRAN"]])
# fortsatt warning og feilmelding,
# Warning in install.packages :
#  'lib = "countrycode"' is not writable
#Error in install.packages : unable to install packages
# kjører på nytt med å svare ja på spørsmålet om
# "personalized library instead", annen feilmelding:
# Warning in install.packages :
#  'lib = "countrycode"' is not writable
#Error in install.packages : argument is not interpretable as logical

# Installsjons fra Package-omr?det virker med f?lgende
# pakkeliste limt inn der:
# Rmpfr needed for crossrun
brew,countrycode,devtools,dplyr,ggplot2,
googleVis,knitr,MCMCpack,repmis,RCurl,
rmarkdown,texreg,tidyr,WDI,
xtable,Zelig,roxygen2,Rmpfr

# laste ned pakkene:
library(brew)
library(countrycode)
library(devtools)
library(dplyr)
library(ggplot2)
library(googleVis)
library(knitr)
library(MCMCpack)
library(repmis)
library(RCurl)
library(rmarkdown)
library(texreg)
library(tidyr)
library(WDI)
library(xtable)
library(Zelig)
library(XML)
library(formatR)
# og roxygen2 og testthat, og pakke for crossrun:
library(roxygen2)
library(testthat)
library(Rmpfr)

# ikke kjøre dette hver gang:
devtools::install_github("r-lib/devtools")
# dette tar lang tid og det hender mye, det kommer også
# opp følgende innimellom:
# * installing *source* package 'git2r' ...
#    **********************************************
# WARNING: this package has a configure script
# It probably needs manual configuration
# **********************************************
# men manual configuration er jo ikke sagt noe sted
# hvordan kan gjøres
# også
# installing to C:/Users/Tore_2/Documents/R/win-library/3.5/git2r/libs/x64
# Warning in file.copy(files, dest, overwrite = TRUE) :
#  problem copying .\git2r.dll to C:\Users\Tore_2\Documents\R\win-library\3.5\git2r\libs\x64\git2r.dll: Permission denied
# som ikke er forwståelig, og
# package ‘stringi’ successfully unpacked and MD5 sums checked
# Warning: cannot remove prior installation of package ‘stringi’
# som heller ikke er mulig å forstå hva kan gjøres med

# opprette pakke:
#devtools::create("./crossrun") # melding om erstattet med:
# usethis::create_package("./crossrun", rstudio=TRUE)
# gir opp å kjøre dette fullt ut reproduserbart ved
# kommandoer i skriptet, men går tilbake til å kjøre
# det fra File, New project og dialogbokser der
# det åpner seg et nytt RStudio.vindu, lukker det
# gamle og åpner lagepakke.r i det nye og fortsetter her

# sjekke om alt som skal være er der:
devtools::has_devel() # skal slutte med TRUE, oks

# kjøre dev_mode fra devtools:
devtools::load_all()
# dev_mode() # da kommer denne:
# Error in read.dcf(con) :
# Line starting 'Run in a series of   ...' is malformed!
# og den er det jeg har lastet ned R fra patch og alt
# annet reinstallert for å bli kvitt. Kommer ikke videre.


# tidy_dir("R") # kommentert vekk etter kjørt en gang

# legge inn dependency i DESCRIPTION:
devtools::use_package(package="Rmpfr")

# kjøre etter å ha lagt inn roxygen2-kommandoer
devtools::document() #nda kommer denne:
# Error in read.dcf(path_desc) :
#  Line starting 'Run in a series of   ...' is malformed!
# kommer ikke videre

# ser på dokumentasjon:
?cumsumm
?cumsummcol
?boxprobt
?exactbin
?clshift
?simclbin
?crossrunsymm
?crossrunbin

# testing (kommentert bort inntil videre):
#use_testthat()
#test()
# det er nok mer innviklet, tar det seinere

# namespace:
search() # viser alle tilgjengelige pakker, og
# .GlobalEnv og Autoloads
# det blir feilmeldinger ved Check hvis ikke alle
# funksjoner eksporteres. Er det en måte rundt det?
# I mellomtida eksporterer jeg alle funksjoner.
# bruker funksjoner fra rmfpr og noen fra stats
# med :: , så tar ikke med import i namespace .

# data:
# kjører ikke dette nå, datafilene er for store:
# henter inn workspace og bruker use_data fra det:
# changes to the article 1 directory and loads data:
setwd("..")
setwd("..")
setwd("./crossrun.art1")
getwd() # ok
load("crossrun1.RData")
# back to package directory:
setwd("..")
setwd("./Pakke/crossrun")
getwd() # ok
# included main data in package:
crossrun100sym <- cr100$pt
use_data(crossrun100sym, overwrite=TRUE)
use_data(crb100..6, overwrite=TRUE)
use_data(crb100..7, overwrite=TRUE)
use_data(crb100..8, overwrite=TRUE)
use_data(crb100..9, overwrite=TRUE)
use_data(cl100simbin, overwrite=TRUE)
# removes other data:
obj1 <- ls()
obj2 <- obj1[!(obj1 %in% c("crossrun100sym","crb100..6","crb100..7","crb100..8",
                     "crb100..9","crb100.2p","cl100simbin"))]
rm(list=obj2)
rm(obj1,obj2)
ls()
# datafilene er store, sjekker komprimring:
tools::checkRdaFiles("/data")
# ikke noe forståelig der, tar vekk data
# og workspace manuelt

build("../crossrun")
# sjekk som cran:

?document

# lage vignett:
devtools::use_vignette("vignetteCrossrun")
# da kommer den forferdelige feilmeldinga
# Error in read.dcf(path_desc) :
# Line starting 'Run in a series of   ...' is malformed!
# igjen. Så dette virker ikke. Lager vignett direkte
# som Rmd-fil.
