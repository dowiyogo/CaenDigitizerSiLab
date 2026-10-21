#ifndef HODO_MAPPING_JSON_H
#define HODO_MAPPING_JSON_H
// Mapeo de canales leído del JSON del DAQ v2 (daqv2.mapping/1), la misma fuente que usa daqrun
// (DAQ2-11): una sola fuente de verdad para las campañas de octubre 2026 en adelante.
//   - Ruta: variable de entorno HODO_MAPPING_JSON o, si no está, el mapeo instalado
//     ~/.local/opt/daqv2/share/mapping/hodo_mapping_2026_10_r3.json.
//   - strip -> canal por plano con la coordenada "strip" del JSON (la medida en r3).
//   - Si el archivo no existe o no se puede leer: kCfg_2026_10 (hoja de octubre) con un aviso en
//     stderr. Las corridas v1-v10 siguen con kCfg_2026_05 (hodo_mapping.h).
// Se incluye desde hodo_mapping.h; necesita nlohmann/json (paquete del sistema).
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>

#include <nlohmann/json.hpp>

namespace hodo_json_detail {
inline std::string default_mapping_path() {
    if (const char* p = std::getenv("HODO_MAPPING_JSON")) return p;
    const char* home = std::getenv("HOME");
    return std::string(home ? home : "") + "/.local/opt/daqv2/share/mapping/hodo_mapping_2026_10_r3.json";
}
}  // namespace hodo_json_detail

// Configuración desde el JSON (se lee una vez por proceso). fallback: la que se devuelve con aviso.
inline const HodoMappingConfig& config_from_mapping_json(const HodoMappingConfig& fallback) {
    static bool loaded = false;
    static bool ok = false;
    static HodoMappingConfig cfg{};
    static std::string name;
    if (!loaded) {
        loaded = true;
        const std::string path = hodo_json_detail::default_mapping_path();
        try {
            std::ifstream in(path);
            if (!in) throw std::runtime_error("no se pudo abrir");
            const auto j = nlohmann::json::parse(in);
            std::array<std::array<int, 16>, 4> t{};
            for (auto& a : t) a.fill(-1);
            for (const auto& e : j.at("channels")) {
                const int ch = e.at("dig_ch").get<int>(), k = e.at("strip").get<int>();
                const bool y = e.at("plane").get<std::string>() == "Y";
                const int st = e.at("station").get<std::string>() == "hod2019" ? 1 : 0;
                if (k >= 0 && k < 16) t[static_cast<std::size_t>(2 * st + (y ? 1 : 0))][static_cast<std::size_t>(k)] = ch;
            }
            for (const auto& a : t)
                for (int ch : a)
                    if (ch < 0) throw std::runtime_error("faltan tiras en un plano");
            name = "json:" + j.at("mapping_version").get<std::string>();
            cfg = HodoMappingConfig{name.c_str(), t[0], t[1], t[2], t[3]};
            ok = true;
        } catch (const std::exception& ex) {
            std::fprintf(stderr, "AVISO hodo_mapping: no se pudo leer el mapeo %s (%s); se usa %s\n", path.c_str(), ex.what(),
                         fallback.name);
        }
    }
    return ok ? cfg : fallback;
}

#endif  // HODO_MAPPING_JSON_H
