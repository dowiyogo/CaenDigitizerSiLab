lib = CaenDigitizerSiLab
SRCS = $(lib).cxx # source files
OBJS = $(SRCS:.cxx=.o)

ROOTCONFIG := root-config
ROOTCFLAGS := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
ROOTGLIBS := $(shell $(ROOTCONFIG) --glibs)

CXX := g++
CXXFLAGS := -O2 -Wall -fPIC $(ROOTCFLAGS)
LD := g++
LDFLAGS := -O2 $(ROOTLDFLAGS)

INCLUDES := -I.
LIBS := $(ROOTGLIBS) -lSpectrum -lEG
HEADERS := $(lib).h

all: lib$(lib).so.1.0.1 link example test_open test_acquire_sw test_acquire_ext_nim reconstruct_v8

reconstruct_v8: reconstruct_v8.cxx
	$(CXX) $(ROOTCFLAGS) -std=c++17 -O2 -Wall reconstruct_v8.cxx \
		$(ROOTLDFLAGS) $(ROOTGLIBS) -o reconstruct_v8

infer_mapping: infer_mapping.cxx
	$(CXX) $(ROOTCFLAGS) -std=c++17 -O2 -Wall infer_mapping.cxx \
		$(ROOTLDFLAGS) $(ROOTGLIBS) -o infer_mapping

reconstruct_v8b: reconstruct_v8b.cxx
	$(CXX) $(ROOTCFLAGS) -std=c++17 -O2 -Wall reconstruct_v8b.cxx \
		$(ROOTLDFLAGS) $(ROOTGLIBS) -o reconstruct_v8b

analyze_highgain: analyze_highgain.cxx
	$(CXX) $(ROOTCFLAGS) -std=c++17 -O2 -Wall analyze_highgain.cxx \
		$(ROOTLDFLAGS) $(ROOTGLIBS) -o analyze_highgain

multirun_analysis: multirun_analysis.cxx
	$(CXX) $(ROOTCFLAGS) -std=c++17 -O2 -Wall multirun_analysis.cxx \
		$(ROOTLDFLAGS) $(ROOTGLIBS) -o multirun_analysis

tagprobe_analysis: tagprobe_analysis.cxx
	$(CXX) $(ROOTCFLAGS) -std=c++17 -O2 -Wall tagprobe_analysis.cxx \
		$(ROOTLDFLAGS) $(ROOTGLIBS) -o tagprobe_analysis

reconstruct_corrected: reconstruct_corrected.cxx
	$(CXX) $(ROOTCFLAGS) -std=c++17 -O2 -Wall reconstruct_corrected.cxx \
		$(ROOTLDFLAGS) $(ROOTGLIBS) -o reconstruct_corrected

reconstruct_cosmic_tracks: reconstruct_cosmic_tracks.cxx
	$(CXX) $(ROOTCFLAGS) -std=c++17 -O2 -Wall reconstruct_cosmic_tracks.cxx \
		$(ROOTLDFLAGS) $(ROOTGLIBS) -o reconstruct_cosmic_tracks

charge_analysis: charge_analysis.cxx
	$(CXX) $(ROOTCFLAGS) -std=c++17 -O2 -Wall charge_analysis.cxx \
		$(ROOTLDFLAGS) $(ROOTGLIBS) -o charge_analysis

poisson_analysis: poisson_analysis.cxx
	$(CXX) $(ROOTCFLAGS) -std=c++17 -O2 -Wall poisson_analysis.cxx \
		$(ROOTLDFLAGS) $(ROOTGLIBS) -o poisson_analysis

analyze_independent: analyze_independent.cxx
	$(CXX) $(ROOTCFLAGS) -std=c++17 -O2 -Wall analyze_independent.cxx \
		$(ROOTLDFLAGS) $(ROOTGLIBS) -lSpectrum -o analyze_independent

example: example.o configuration.o CaenDigitizerSiLab.o CaenDigitizerSiLabDict.o
	$(CXX) $(LDFLAGS) $^ -L. -lCaenDigitizerSiLab -lCAENDigitizer $(LIBS) -o $@

test_open: test_open.o configuration.o CaenDigitizerSiLab.o CaenDigitizerSiLabDict.o
	$(CXX) $(LDFLAGS) $^ -L. -lCaenDigitizerSiLab -lCAENDigitizer $(LIBS) -o $@

test_acquire_sw: test_acquire_sw.o configuration.o CaenDigitizerSiLab.o CaenDigitizerSiLabDict.o
	$(CXX) $(LDFLAGS) $^ -L. -lCaenDigitizerSiLab -lCAENDigitizer $(LIBS) -o $@

test_acquire_ext_nim: test_acquire_ext_nim.o configuration.o CaenDigitizerSiLab.o CaenDigitizerSiLabDict.o
	$(CXX) $(LDFLAGS) $^ -L. -lCaenDigitizerSiLab -lCAENDigitizer $(LIBS) -o $@

example.o : example.cxx
	$(CXX) $(INCLUDES) $(CXXFLAGS) -c $^ -o $@

test_open.o : test_open.cxx
	$(CXX) $(INCLUDES) $(CXXFLAGS) -c $^ -o $@

test_acquire_sw.o : test_acquire_sw.cxx
	$(CXX) $(INCLUDES) $(CXXFLAGS) -c $^ -o $@

test_acquire_ext_nim.o : test_acquire_ext_nim.cxx
	$(CXX) $(INCLUDES) $(CXXFLAGS) -c $^ -o $@

configuration.o: configuration.C
	$(CXX) $(INCLUDES) $(CXXFLAGS) -c $^ -o $@
	#$(CXX) -c -Wall -DLINUX -DUNIX $^ -o $@
	#$(CXX) $(INCLUDES) -Wall -DLINUX -DUNIX -c $^ -o $@

$(lib)Dict.cxx: configuration.h $(HEADERS) LinkDef.h
	rootcling -f $(lib)Dict.cxx $^

CaenDigitizerSiLabDict.o: CaenDigitizerSiLabDict.cxx
	$(CXX) $(CXXFLAGS) -c $^ -o $@

lib$(lib).so.1.0.1 : $(OBJS) CaenDigitizerSiLabDict.o
	$(CXX) $(LDFLAGS) $(LIBS) $(CXXFLAGS) -shared -Wl,-soname,lib$(lib).so.1 $^ -lc -o $@

%.o : %.cxx
	$(CXX) $(CXXFLAGS) -g -c -Wall -fPIC $^ -o $@

.PHONY: CLEAN link
clean:
	rm -f *.o *.d *.pcm *.root *.txt CaenDigitizerSiLabDict.cxx CaenDigitizerSiLabDict.o example test_open test_acquire_sw test_acquire_ext_nim

link:
	ln -sf libCaenDigitizerSiLab.so.1.0.1 libCaenDigitizerSiLab.so.1
	ln -sf libCaenDigitizerSiLab.so.1 libCaenDigitizerSiLab.so
