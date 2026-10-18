// analysis/src/hodo_paths.h
// Centralized path management for raw data in the CAEN digitizer analysis pipeline.
// Conforms to TASK PATHS-1 specification.

#ifndef HODO_PATHS_H
#define HODO_PATHS_H

#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

namespace hodo {

// Returns the raw data directory path without a trailing slash.
// Uses environment variable HODO_RAW_DIR if set and non-empty;
// otherwise defaults to "/home/xplus/digitizer/data/raw".
inline std::string rawDataDir() {
    const char* env = std::getenv("HODO_RAW_DIR");
    if (env != nullptr && env[0] != '\0') {
        std::string s(env);
        while (s.size() > 1 && s.back() == '/') {
            s.pop_back();
        }
        return s;
    }
    return "/home/xplus/digitizer/data/raw";
}

// Resolves a raw data filename or path:
// 1. If HODO_RAW_DIR is set and non-empty, prioritize looking inside rawDataDir().
//    Collects all matches (direct and in each first-level subdirectory). If more
//    than one match is found, prints an error to stderr listing all matches and returns "".
// 2. Otherwise, check if fname exists as given (e.g. local symlink or relative/absolute path).
// 3. Otherwise, check default rawDataDir() and its direct subdirectories. If more
//    than one match is found, prints an error to stderr listing all matches and returns "".
// 4. Return the resolved path if found; otherwise return fname as fallback.
inline std::string resolveRawFile(const std::string& fname) {
    if (fname.empty()) return fname;

    namespace fs = std::filesystem;
    const char* env = std::getenv("HODO_RAW_DIR");
    if (env != nullptr && env[0] != '\0') {
        fs::path base(rawDataDir());
        std::vector<std::string> matches;
        std::error_code ec;
        fs::path p_direct = base / fname;
        if (fs::exists(p_direct, ec)) {
            matches.push_back(p_direct.string());
        }
        if (fs::is_directory(base, ec)) {
            for (const auto& entry : fs::directory_iterator(base, ec)) {
                if (entry.is_directory(ec)) {
                    fs::path p_sub = entry.path() / fname;
                    if (fs::exists(p_sub, ec)) {
                        matches.push_back(p_sub.string());
                    }
                }
            }
        }
        if (matches.size() > 1) {
            std::cerr << "ERROR: Ambiguous raw file match for '" << fname << "' in base directory " << base.string() << ":\n";
            for (const auto& m : matches) {
                std::cerr << "  " << m << "\n";
            }
            return "";
        }
        if (matches.size() == 1) {
            return matches[0];
        }
    }

    std::error_code ec;
    if (fs::exists(fname, ec)) {
        return fname;
    }

    fs::path base_default(rawDataDir());
    std::vector<std::string> matches;
    fs::path p_def_direct = base_default / fname;
    if (fs::exists(p_def_direct, ec)) {
        matches.push_back(p_def_direct.string());
    }
    if (fs::is_directory(base_default, ec)) {
        for (const auto& entry : fs::directory_iterator(base_default, ec)) {
            if (entry.is_directory(ec)) {
                fs::path p_sub = entry.path() / fname;
                if (fs::exists(p_sub, ec)) {
                    matches.push_back(p_sub.string());
                }
            }
        }
    }
    if (matches.size() > 1) {
        std::cerr << "ERROR: Ambiguous raw file match for '" << fname << "' in base directory " << base_default.string() << ":\n";
        for (const auto& m : matches) {
            std::cerr << "  " << m << "\n";
        }
        return "";
    }
    if (matches.size() == 1) {
        return matches[0];
    }

    return fname;
}

} // namespace hodo

// Provide rawDataDir in global namespace for convenience
inline std::string rawDataDir() {
    return hodo::rawDataDir();
}

#endif // HODO_PATHS_H

