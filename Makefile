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

all: lib$(lib).so.1.0.1 link example test_open test_acquire_sw test_acquire_ext_nim

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
	rm -f *.o *.d *.pcm CaenDigitizerSiLabDict.cxx CaenDigitizerSiLabDict.o example test_open test_acquire_sw test_acquire_ext_nim

link:
	ln -sf libCaenDigitizerSiLab.so.1.0.1 libCaenDigitizerSiLab.so.1
	ln -sf libCaenDigitizerSiLab.so.1 libCaenDigitizerSiLab.so
