#############################################################
###  User definitions section:

###  Postprocessing calibration analysis program names
CALIBRATIONANALYZERS := ResidualsCalculator
###  Generated input XML files to Shine
GENERATEDXMLS := $(patsubst %.xml.in,%.xml,$(wildcard *.xml.in))

###  Symbolic targets (not real generated files)
.PHONY: tools clean

###  Generate necessary input XMLs and programs
tools: $(GENERATEDXMLS) $(CALIBRATIONANALYZERS)
	@chmod +x RunResidualsCalculator.sh

#############################################################
###  The parts below are for experts only (steering the compilation):
###  set up input XML configuration and compile analyzer programs

SHINEOFFLINECONFIG := shine-offline-config
CONFIGFILES := $(shell $(SHINEOFFLINECONFIG) --config)
DBFILES := $(shell $(SHINEOFFLINECONFIG) --db-path)
DOCPATH := $(shell $(SHINEOFFLINECONFIG) --doc)
CFLAGS := $(shell $(SHINEOFFLINECONFIG) --cppflags)
CPPFLAGS := $(shell $(SHINEOFFLINECONFIG) --cppflags)
CXXFLAGS := $(shell $(SHINEOFFLINECONFIG) --cxxflags)
LDFLAGS  := $(shell $(SHINEOFFLINECONFIG) --ldflags)

ifeq "$(strip $(CXX))" ""
  CXX := g++
endif

###  For generating input XMLs
%: %.in
	@echo -n "Generating $@ file..."
	@ sed -e 's!@''CONFIGDIR@!$(CONFIGFILES)!g;s!@''SHINEDBDIR@!$(DBFILES)!g' $< >$@
	@echo "done"
	@chmod +x RunResidualsCalculator.sh

###  For compiling analyzer programs
%: %.C
	$(CXX) -o $@ $< `root-config --cflags` `root-config --libs` $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS)
%: %.cpp
	$(CXX) -o $@ $< `root-config --cflags` `root-config --libs` $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS)
%: %.cc
	$(CXX) -o $@ $< `root-config --cflags` `root-config --libs` $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS)
%: %.cxx
	$(CXX) -o $@ $< `root-config --cflags` `root-config --libs` $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS)

###  For cleaning up
clean:
	rm -f $(GENERATEDXMLS) $(CALIBRATIONANALYZERS) $(CALIBRATIONTARGETS) *.root *.pdf produce merge calibrate run
