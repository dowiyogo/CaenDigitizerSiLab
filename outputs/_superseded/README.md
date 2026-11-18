# Superseded Analysis Outputs

This directory stores historical analysis outputs that have been superseded by runs using the verified authoritative physical detector mapping.

- `track_reconstruction_interleaved`: Superseded by `track_reconstruction_interleaved_v7_corrected_mapping` because its `summary.txt` documents an uncorrected sequential channel mapping (hod2019 X0=Ch48, Y0=Ch32) that differs from the authoritative physical cabling `kH19X[0]=53, kH19Y[0]=37` [V: `summary.txt:18-19`, `reconstruct_corrected.cxx:62-63`].
- `v7_track_reconstruction`: Superseded by `track_reconstruction_interleaved_v7_corrected_mapping` because its `summary.txt` documents an uncorrected sequential channel mapping (hod2019 X0=Ch48, Y0=Ch32) that differs from the authoritative physical cabling `kH19X[0]=53, kH19Y[0]=37` [V: `summary.txt:18-19`, `reconstruct_corrected.cxx:62-63`].

- `2026-11-04_followup14/oct_boot_fit_0_0_empty.csv`: archivo vacío de la compilación de boot_fit.C (filas 0 a 0), sin datos.
- `2026-11-19_response_matrix_summary_v1/`: primera salida de `rm_summary.py` (matriz de respuesta, 2026-11-19). La comprobación de 1 p.e. usaba una gaussiana sola para el pico oscuro de n = 0, que no converge sobre el continuo de ruido (posiciones negativas); las matrices son las mismas que en la versión vigente de `outputs/hod2018_response_matrix/2026-11-19/`.
- `2026-11-19_response_matrix_core_v1/`: primera salida de `rm_core.py` (2026-11-19), solo con las variantes ls y gnorm y con un núcleo sin sentido en n = 60 cuando el desborde corta la ventana del ajuste; reemplazada por la versión con las cuatro variantes (`window_edge.csv` se regenera igual).
