sources := $(wildcard *.tex)
plots := $(wildcard */*.eps)

%.dvi: %.tex $(plots)
	latex $< ; latex $<

%.ps: %.dvi
	dvips $< -o $@

%.pdf: %.dvi
	dvipdf $< $@

note : $(sources:.tex=.ps)

pdf : $(sources:.tex=.pdf)

.PHONY: all clean edit view viewpdf 

all : note

clean:
	rm -f *.dvi *.ps *.aux *.toc *.log *.pdf *~

view: note
	gv ztomumu.ps &

viewpdf: pdf
	xpdf ztomumu.pdf &

edit:
	emacs ztomumu.tex &

