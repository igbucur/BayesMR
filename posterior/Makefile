POLYCHORD_DIR = ../inst/PolyChordLite
LIB_DIR = $(POLYCHORD_DIR)/lib
DRIVER = polychord_BayesMR
LIKELIHOOD = BayesMR_likelihood
LDLIBS = -lstdc++ -larmadillo

default: $(DRIVER).o lib$(LIKELIHOOD).a
	$(FC) $(DRIVER).o -o ../BayesMR -L. -L$(POLYCHORD_DIR)/lib -Wl,-Bstatic -lchord -Wl,-Bdynamic -l$(LIKELIHOOD) $(LDLIBS)

# rule for building the PolyChord driver
$(DRIVER).o: $(DRIVER).cpp

# rule for building the c++ wrapper
lib$(LIKELIHOOD).a: $(LIKELIHOOD).o
	$(AR) $@ $^ 

# Utility targets
.PHONY: clean veryclean

clean:
	$(RM) *.o *.mod *.MOD

veryclean: clean
	$(RM) *~ *.a
