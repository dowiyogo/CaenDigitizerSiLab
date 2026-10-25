/* test_sy5527.c
 * Prueba mínima de conexión al SY5527 vía CAEN HV Wrapper (TCP/IP).
 * Uso: ./test_sy5527 <IP> [usuario] [password]
 *      usuario y password salen de los argumentos o, si faltan, de DAQ_HV_USER y DAQ_HV_PASSWORD
 *      en el entorno (como hvctl); sin valores por defecto. Reemplazado por `hvctl map`.
 * Salida: código de retorno de la Wrapper y mapa de slots del crate.
 * Compilar: gcc -Wall test_sy5527.c -o test_sy5527 -lcaenhvwrapper
 * Nota: verificar los nombres SY5527 y LINKTYPE_TCPIP en el CAENHVWrapper.h instalado.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "CAENHVWrapper.h"

int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "Uso: %s <IP> [usuario] [password]\n", argv[0]);
        return 1;
    }
    const char *ip   = argv[1];
    const char *user = (argc > 2) ? argv[2] : getenv("DAQ_HV_USER");
    const char *pass = (argc > 3) ? argv[3] : getenv("DAQ_HV_PASSWORD");
    if (!user || !*user || !pass || !*pass) {
        fprintf(stderr, "Faltan usuario o password: argumentos o DAQ_HV_USER/DAQ_HV_PASSWORD\n");
        return 1;
    }
    int handle = -1;

    /* SY5527: tipo de sistema; con LINKTYPE_TCPIP el argumento es la IP como string */
    CAENHVRESULT ret = CAENHV_InitSystem(SY5527, LINKTYPE_TCPIP,
                                         (void *)ip, user, pass, &handle);
    if (ret != CAENHV_OK) {
        fprintf(stderr, "InitSystem fallo: codigo 0x%X -> %s\n",
                (unsigned)ret, CAENHV_GetError(handle));
        return 1;
    }
    printf("Conectado a %s. handle = %d\n", ip, handle);

    unsigned short nSlots, *nCh, *serial;
    unsigned char *fwMin, *fwMax;
    char *models, *descr;

    ret = CAENHV_GetCrateMap(handle, &nSlots, &nCh, &models, &descr,
                             &serial, &fwMin, &fwMax);
    if (ret == CAENHV_OK) {
        char *m = models, *d = descr;   /* strings concatenados separados por '\0' */
        for (int s = 0; s < nSlots; s++) {
            if (*m)
                printf("Slot %2d: %-8s %-30s ch=%u SN=%u FW=%u.%02u\n",
                       s, m, d, nCh[s], serial[s], fwMax[s], fwMin[s]);
            else
                printf("Slot %2d: vacio\n", s);
            m += strlen(m) + 1;
            d += strlen(d) + 1;
        }
        CAENHV_Free(nCh);   CAENHV_Free(serial);
        CAENHV_Free(fwMin); CAENHV_Free(fwMax);
        CAENHV_Free(models); CAENHV_Free(descr);
    } else {
        fprintf(stderr, "GetCrateMap fallo: %s\n", CAENHV_GetError(handle));
    }

    CAENHV_DeinitSystem(handle);
    return 0;
}
