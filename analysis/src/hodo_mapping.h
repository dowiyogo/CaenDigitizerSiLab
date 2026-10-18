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

// Strip k -> channel mappings
static constexpr std::array<int, 16> kH18X = {15, 0, 13, 2, 11, 4, 9, 6, 7, 8, 5, 10, 3, 12, 1, 14};
static constexpr std::array<int, 16> kH18Y = {16, 24, 17, 25, 18, 26, 19, 27, 20, 28, 21, 29, 22, 30, 23, 31};
static constexpr std::array<int, 16> kH19X = {53, 61, 54, 62, 55, 63, 52, 60, 51, 59, 50, 58, 49, 57, 48, 56};
static constexpr std::array<int, 16> kH19Y = {37, 45, 38, 46, 39, 47, 36, 44, 35, 43, 34, 42, 33, 41, 32, 40};

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
