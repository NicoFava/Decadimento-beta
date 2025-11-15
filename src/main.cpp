#include "funzioni.h"

int main() {
    TApplication app("app", nullptr, nullptr); // Inizializza l'applicazione ROOT

    run_dalitz();

    app.Run(); // Mostra i canvas e mantiene l'applicazione in esecuzione

    return 0;
}