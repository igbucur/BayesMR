
default: build_PolyChord polychord_MR

build_PolyChord: 
	cd polychordlite && $(MAKE) pypolychord && python setup.py -q install --user

polychord_MR:
	cd polychordlite && $(MAKE) polychord_MR
	cp polychordlite/bin/polychord_MR bin/polychord_MR

clean_MR:
	cd polychordlite/likelihoods/MR && rm *.o
	cd polychordlite/lib && rm libMR_likelihood.a
	cd polychordlite/src/drivers && rm polychord_MR.o

clean: 
	cd polychordlite && $(MAKE) clean

veryclean:
	cd polychordlite && $(MAKE) veryclean



