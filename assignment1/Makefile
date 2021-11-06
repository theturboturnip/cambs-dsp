.PHONY: all html-pdf
all: html-pdf

clean:
	mkdir -p ./build/
	rm -r ./build/*

./build/ex_complete.m: ex_1a.m ex_1b.m ex_1c.m ex_1d.m covid_new_cases_2021_11_04.csv
	mkdir -p ./build/
	cat ex_1a.m ex_1b.m ex_1c.m ex_1d.m > ./build/ex_complete.m
	cp covid_new_cases_2021_11_04.csv ./build/

./build/html/ex_complete.tex: ./build/ex_complete.m
	@echo "Building LaTeX from MATLAB..."
	(cd ./build/ && matlab -batch "publish('./ex_complete.m', 'latex'); exit")
	@echo "...done!"

./build/html/ex_complete.html: ./build/ex_complete.m
	@echo "Building HTML from MATLAB..."
	(cd ./build/ && matlab -batch "publish('./ex_complete.m', 'html'); exit")
	@echo "...done!"

matlab-pdf: ./ex_complete_matlab.pdf
./ex_complete_matlab.pdf: ./build/ex_complete.m
	@echo "Building PDF from MATLAB..."
	(cd ./build/ && matlab -batch "publish('./ex_complete.m', 'pdf'); exit")
	@echo "...done!"
	cp ./build/html/ex_complete.pdf ./ex_complete_matlab.pdf

latex-pdf: ./ex_complete_latex.pdf
./ex_complete_latex.pdf: ./build/html/ex_complete.tex
	(cd ./build/html/ && pdflatex ./ex_complete.tex)
	cp ./build/html/ex_complete.pdf ./ex_complete_latex.pdf

html-pdf: ./build/html/ex_complete.html
	(cd ./build/html/ && google-chrome-stable --new-window --print-to-pdf ./ex_complete.html)
