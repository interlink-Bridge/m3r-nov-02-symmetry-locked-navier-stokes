PDF=paper/main.pdf
TEX=paper/main.tex

pdf:
	cd paper && latexmk -pdf main.tex

clean:
	cd paper && latexmk -C
