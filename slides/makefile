RMDFILES = $(wildcard *.Rmd)

all: slides.html fdat_cpr_run.html

slides.html: $(RMDFILES) style.css template.html render-slides.R
	R CMD BATCH --vanilla render-slides.R

fdat_cpr_run.html: fdat_cpr_run.R
	R CMD BATCH --vanilla fdat_cpr_run.R

