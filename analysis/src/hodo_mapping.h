#ifndef HODO_MAPPING_H
#define HODO_MAPPING_H

#include <array>

// =============================================================================
// Hodoscope Channel-to-Strip Mapping
// =============================================================================
//
// Sources:
// - hod2018X and hod2018Y: Derived in TASK: ANALYSIS-FIX-8 (2026-10-18) from
//   spatial neighbor coincidences (excess of same-event over mixed-event pairs
//   measured in HOD-DIAG-1, outputs/hod_channel_correlations_20261017_202452).
//   The spatial ordering makes |k_a - k_b| = 1 for every top coincident pair.
//   The overall orientation (strip k vs 15-k) is open (not uniquely determined).
//   The position of dead channel 3 is inferred by elimination and continuity.
//
// - hod2019X and hod2019Y: Authoritative cabling record "Disposición Correcta"
//   [V: reconstruct_corrected.cxx:62-63]. Results produced with sequential
//   mapping (X0=48, Y0=32) are superseded.
//
// Coordinate conventions:
// - Strip index k runs 0..15 for each of the 4 planes.
// - Physical pitch = 2.0625 mm (33.0 mm / 16 bars) or 2.0 mm nominal.
// =============================================================================

// Planes: 0=18X, 1=18Y, 2=19X, 3=19Y
constexpr int kPlane18X = 0;
constexpr int kPlane18Y = 1;
constexpr int kPlane19X = 2;
constexpr int kPlane19Y = 3;
// Configuration container
struct HodoMappingConfig {
    const char* name;
    std::array<int, 16> h18x;
    std::array<int, 16> h18y;
    std::array<int, 16> h19x;
    std::array<int, 16> h19y;
};

// Configuration kCfg_2026_05:
// hod2018 from neighbour correlations, hod2019 from sheet "Mayo 2026".
// Note: the sheet "Mayo 2026" gives hod2018 as ch k, which the neighbour correlations of run v8 contradict.
static constexpr HodoMappingConfig kCfg_2026_05 = {
    "kCfg_2026_05",
    {15, 0, 13, 2, 11, 4, 9, 6, 7, 8, 5, 10, 3, 12, 1, 14},
    {16, 24, 17, 25, 18, 26, 19, 27, 20, 28, 21, 29, 22, 30, 23, 31},
    {53, 61, 54, 62, 55, 63, 52, 60, 51, 59, 50, 58, 49, 57, 48, 56},
    {37, 45, 38, 46, 39, 47, 36, 44, 35, 43, 34, 42, 33, 41, 32, 40}
};
// HV-trip check (runs v8, v10): hod2019 faces -> digitizer groups: Delantera 48-55, Trasera 56-63,
// Frontal 40-47, Opuesta 32-39; within a face, connector channel 1..8 -> offset {0,1,2,3,4,7,6,5}.
// 15 of 16 predicted dead channels observed (Opuesta 8 -> ch37 not observed dead).
// hod2018 from HV trips: ch6 = X Canal 7, ch0 = X Canal 13, ch20 = Y Canal 7, {ch22, ch31} = Y Canal {11, 15},
// {ch1, ch3} = X Canal {0, 2}. Inconsistent with the sheet "Mayo 2026" (ch k = coordinate k) for hod2018.

// Configuration kCfg_2026_10:
// Cabling of "Disposicion Correcta Hodoscopios", sheet "Octubre 2026":
//   hod2018 X strip k -> ch k (k=0..15), hod2018 Y strip k -> ch 16+k,
//   hod2019 X strip k -> ch 32+k, hod2019 Y strip k -> ch 48+k.
static constexpr HodoMappingConfig kCfg_2026_10 = {
    "kCfg_2026_10",
    {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15},
    {16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31},
    {32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47},
    {48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63}
};

// Backward-compatibility strip k -> channel mappings (programs keep using kCfg_2026_05)
static constexpr auto kH18X = kCfg_2026_05.h18x;
static constexpr auto kH18Y = kCfg_2026_05.h18y;
static constexpr auto kH19X = kCfg_2026_05.h19x;
static constexpr auto kH19Y = kCfg_2026_05.h19y;

// Returns configuration for a given run (v1..v10 return kCfg_2026_05)
inline const HodoMappingConfig& config_for_run(int /*run*/) {
    return kCfg_2026_05;
}

inline const HodoMappingConfig& config_for_run(const std::string& /*run*/) {
    return kCfg_2026_05;
}

namespace detail {
    inline std::array<int, 64> make_inv_map(const std::array<int, 16>& arr) {
        std::array<int, 64> inv{};
        inv.fill(-1);
        for (int k = 0; k < 16; ++k) {
            inv[arr[k]] = k;
        }
        return inv;
    }
}

// Inverse mapping: channel -> strip k (or -1 if not in that plane)
inline const std::array<int, 64>& get_inv18x() {
    static const auto inv = detail::make_inv_map(kH18X);
    return inv;
}
inline const std::array<int, 64>& get_inv18y() {
    static const auto inv = detail::make_inv_map(kH18Y);
    return inv;
}
inline const std::array<int, 64>& get_inv19x() {
    static const auto inv = detail::make_inv_map(kH19X);
    return inv;
}
inline const std::array<int, 64>& get_inv19y() {
    static const auto inv = detail::make_inv_map(kH19Y);
    return inv;
}

// Return plane index 0..3 for channel c (0=18X, 1=18Y, 2=19X, 3=19Y, -1 if invalid)
inline int ch_to_plane(int c) {
    if (c >= 0 && c <= 15) return kPlane18X;
    if (c >= 16 && c <= 31) return kPlane18Y;
    if (c >= 48 && c <= 63) return kPlane19X;
    if (c >= 32 && c <= 47) return kPlane19Y;
    return -1;
}

// Return physical strip index k (0..15) for channel c (-1 if invalid)
inline int ch_to_strip(int c) {
    if (c >= 0 && c <= 15) return get_inv18x()[c];
    if (c >= 16 && c <= 31) return get_inv18y()[c];
    if (c >= 48 && c <= 63) return get_inv19x()[c];
    if (c >= 32 && c <= 47) return get_inv19y()[c];
    return -1;
}

// Compatibility helper functions
inline int hod2018_x_index(int ch) { return (ch >= 0 && ch <= 15) ? get_inv18x()[ch] : -1; }
inline int hod2018_y_index(int ch) { return (ch >= 16 && ch <= 31) ? get_inv18y()[ch] : -1; }
inline int hod2019_x_index(int ch) { return (ch >= 48 && ch <= 63) ? get_inv19x()[ch] : -1; }
inline int hod2019_y_index(int ch) { return (ch >= 32 && ch <= 47) ? get_inv19y()[ch] : -1; }

#endif // HODO_MAPPING_H

