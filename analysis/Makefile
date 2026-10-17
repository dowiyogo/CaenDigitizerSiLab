CXX := g++
CXXFLAGS := -O2 -Wall -std=c++17 $(shell root-config --cflags)
LDFLAGS := $(shell root-config --ldflags)
LIBS := $(shell root-config --glibs) -lSpectrum

BIN_DIR := bin
SRC_DIR := src

PROGRAMS := \
	analyze_highgain \
	analyze_hod2018 \
	analyze_independent \
	channel_health_analysis \
	channel_health_refit \
	charge_analysis \
	compare_hv_runs \
	final_hodoscope_analysis \
	infer_mapping \
	multirun_analysis \
	plot_hod2018_results \
	poisson_analysis \
	reconstruct_corrected \
	reconstruct_cosmic_tracks \
	reconstruct_tracks_interleaved \
	reconstruct_v8 \
	reconstruct_v8b \
	scan_active_channels \
	scan_active_channels_windowed \
	spe_gain_analysis \
	tagprobe_analysis

TARGETS := $(addprefix $(BIN_DIR)/, $(PROGRAMS))

.PHONY: all clean

all: $(TARGETS)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

$(BIN_DIR)/%: $(SRC_DIR)/%.cpp | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< $(LDFLAGS) $(LIBS) -o $@

$(BIN_DIR)/%: $(SRC_DIR)/%.cxx | $(BIN_DIR)
	$(CXX) $(CXXFLAGS) $< $(LDFLAGS) $(LIBS) -o $@

clean:
	rm -rf $(BIN_DIR)/*

