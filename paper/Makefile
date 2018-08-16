SRC = parareal

all:
	lualatex $(SRC)
	bibtex $(SRC)
	lualatex $(SRC)
	lualatex $(SRC)

clean:
	rm *.pdf $(SRC).log $(SRC).out $(SRC).aux $(SRC).bbl $(SRC).blg 
