
default: polychord_MR

polychord_MR:
	mkdir -p bin
	mkdir -p chains/clusters
	mkdir -p ini
	cd inst/PolyChordLite && $(MAKE) polychord_MR MPI=
	cp inst/PolyChordLite/bin/polychord_MR bin/polychord_MR

# This only cleans the MR application
clean:
	cd inst/PolyChordLite/likelihoods/MR && rm -rf *.o
	cd inst/PolyChordLite/lib && rm libMR_likelihood.a
	cd inst/PolyChordLite/src/drivers && rm polychord_MR.o

veryclean:
	cd inst/PolyChordLite && $(MAKE) veryclean



