# Genera un encabezado C++ con el contenido de un archivo de texto como cadena cruda.
#   cmake -DIN=<archivo> -DOUT=<encabezado> -DNAME=<identificador> -P embed_text.cmake
file(READ "${IN}" content)
string(FIND "${content}" ")DAQEMBED\"" bad)
if(NOT bad EQUAL -1)
  message(FATAL_ERROR "${IN} contiene el delimitador )DAQEMBED\"")
endif()
file(WRITE "${OUT}.tmp" "// Generado desde ${IN} por embed_text.cmake: no editar.\n#pragma once\ninline constexpr const char* ${NAME} = R\"DAQEMBED(${content})DAQEMBED\";\n")
file(COPY_FILE "${OUT}.tmp" "${OUT}" ONLY_IF_DIFFERENT)
file(REMOVE "${OUT}.tmp")
