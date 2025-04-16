lib = CaenDigitizerSiLab
SRCS = $(lib).cxx
OBJS = $(SRCS:.cxx=.o)

ROOTCONFIG := root-config
ROOTCFLAGS := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
ROOTLIBS := $(shell $(ROOTCONFIG) --libs)

CXX := g++
CXXFLAGS := -O2 -Wall -fPIC $(ROOTCFLAGS)
LD := g++
LDFLAGS := -O2 $(ROOTLDFLAGS)

INCLUDES := -I.
LIBS := -lCore -lRIO $(ROOTLIBS) -lMathMore -lSpectrum -lEG -lCAENDigitizer
HEADERS := $(lib).h

SO_NAME := lib$(lib).so.1.0.1
SO_LINK1 := lib$(lib).so.1
SO_LINK2 := lib$(lib).so

all: $(SO_NAME) install-links

example: $(SO_NAME) example.o
	$(CXX) $(LDFLAGS) -L. example.o -l$(lib) $(LIBS) -o $@

example.o : example.cxx
	$(CXX) $(INCLUDES) $(CXXFLAGS) -c $< -o $@

$(lib)Dict.cxx: $(HEADERS) LinkDef.h
	rootcling -f $(lib)Dict.cxx -s $(SO_NAME) -c $(HEADERS) LinkDef.h

$(SO_NAME): $(OBJS) $(lib)Dict.cxx
	$(CXX) $(LDFLAGS) $(CXXFLAGS) -shared -Wl,-soname,$(SO_LINK1) $(OBJS) $(lib)Dict.cxx $(LIBS) -lc -o $@

%.o : %.cxx
	$(CXX) $(CXXFLAGS) -g -c -Wall -fPIC $< -o $@

install-links:
	ln -sf $(SO_NAME) $(SO_LINK1)
	ln -sf $(SO_LINK1) $(SO_LINK2)

install: all
	sudo cp $(SO_NAME) $(SO_LINK1) $(SO_LINK2) /usr/local/lib/
	sudo ldconfig
	@echo ">> Installed to /usr/local/lib and updated ldconfig."

clean:
	rm -f *.o *.so* $(lib)Dict.* example

