
default: build_PolyChord polychord_MR

build_PolyChord: 
	cd polychordlite && $(MAKE) pypolychord && python setup.py -q install --user

polychord_MR:
	cd polychordlite && $(MAKE) polychord_MR
	cp polychordlite/bin/polychord_MR bin/polychord_MR

# This only cleans the MR application, but it does it properly, unlike PolyChord
clean:
	cd polychordlite/likelihoods/MR && rm -rf *.o
	cd polychordlite/lib && rm libMR_likelihood.a
	cd polychordlite/src/drivers && rm polychord_MR.o

veryclean:
	cd polychordlite && $(MAKE) veryclean



