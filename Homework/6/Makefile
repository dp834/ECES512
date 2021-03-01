MATHEMATICA_FILES := $(wildcard *.wls)
BUILD_WHEN_NEEDED := $(addprefix .mod_,$(MATHEMATICA_FILES))

MATLAB_FILES = Homework_6.m
MATLAB_MAIN  = Homework_6.m

IMG_OUTPUT = images/

LATEX_SRC  = Homework-6.tex
PDF_OUTPUT = Homework-6

run: $(MATLAB_FILES) $(MATLAB_MAIN)
	mkdir -p $(IMG_OUTPUT)
	matlab -nodisplay -nosplash -batch "run('$(MATLAB_MAIN)')"

math: $(BUILD_WHEN_NEEDED)

$(BUILD_WHEN_NEEDED): .mod_% : %
	@echo "Running $^"
	@math -noprompt -run "<<$^; Exit[]"
	@touch $@

pdf: $(LATEX_SRC)
	pdflatex -jobname $(PDF_OUTPUT) $(LATEX_SRC)

clean:
	@rm -rf .mod_* $(IMG_OUTPUT)/* $(PDF_OUTPUT).pdf *.aux *.log *.nav *.toc *.snm *.out &2> /dev/null || true
