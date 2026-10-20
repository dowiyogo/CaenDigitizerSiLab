// dash_server.h
// THttpServer que además sirve una página estática de solo números en /dash/ (también /dash y
// /dash/index.html). El HTML va embebido en el binario (generado en la compilación desde
// apps/<app>/dash.html), así que no depende de archivos instalados. El resto de las peticiones
// (gráficos, /status/...) sigue en THttpServer sin cambios. Los archivos estáticos de civetweb no
// sirven aquí: con AddLocation, "/dash/" da 403 (listado de directorio).
// Seguridad entre hilos: como THttpServer (ProcessRequest corre en el hilo que llama ProcessRequests).
#pragma once

#include <memory>
#include <string>
#include <utility>

#include <THttpCallArg.h>
#include <THttpServer.h>

namespace daq::io {

class DashServer : public THttpServer {
 public:
  DashServer(const char* engine, std::string html) : THttpServer(engine), html_(std::move(html)) {}
  // Ruta de la petición (sin la barra inicial) que corresponde a la página.
  static bool isDashRequest(const std::string& path, const std::string& file) {
    return ((path == "dash" || path == "dash/") && (file.empty() || file == "index.html")) || (path.empty() && file == "dash");
  }

 protected:
  void ProcessRequest(std::shared_ptr<THttpCallArg> arg) override {
    if (isDashRequest(arg->GetPathName(), arg->GetFileName())) {
      arg->SetContentType("text/html; charset=utf-8");
      arg->AddNoCacheHeader();
      arg->SetContent(std::string(html_));
      return;
    }
    THttpServer::ProcessRequest(arg);
  }

 private:
  std::string html_;
};

}  // namespace daq::io
