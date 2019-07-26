
default: build_PolyChord polychord_MR

build_PolyChord: 
	cd inst/PolyChordLite && $(MAKE) pypolychord && python setup.py -q install --user

polychord_MR:
	mkdir -p bin
	mkdir -p chains/clusters
	mkdir -p ini
	cd inst/PolyChordLite && $(MAKE) polychord_MR
	cp inst/PolyChordLite/bin/polychord_MR bin/polychord_MR

# This only cleans the MR application, but it does it properly, unlike PolyChord
clean:
	cd inst/PolyChordLite/likelihoods/MR && rm -rf *.o
	cd inst/PolyChordLite/lib && rm libMR_likelihood.a
	cd inst/PolyChordLite/src/drivers && rm polychord_MR.o

veryclean:
	cd inst/PolyChordLite && $(MAKE) veryclean



