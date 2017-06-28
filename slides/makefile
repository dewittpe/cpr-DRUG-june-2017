RMDFILES = $(wildcard *.Rmd)

all: slides.html fdat_cpr_run.html influence-examples.Rout

slides.html: $(RMDFILES) style.css template.html render-slides.R
	R CMD BATCH --vanilla render-slides.R

fdat_cpr_run.html: fdat_cpr_run.R
	R CMD BATCH --vanilla fdat_cpr_run.R

influence-examples.Rout: influence-examples.R
	R CMD BATCH --vanilla influence-examples.R
