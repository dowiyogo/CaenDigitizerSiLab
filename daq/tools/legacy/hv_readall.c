/* hv_readall.c
 * Lectura (SOLO LECTURA, ninguna llamada Set) de todos los parámetros de canal
 * de una tarjeta del SY5527 vía CAEN HV Wrapper.
 *
 * Uso:    ./hv_readall <IP> <slot> [usuario] [password]
 *         ./hv_readall 192.168.1.118 3 admin admin
 * Salida: lista de parámetros del canal 0 con su tipo y unidad, y una tabla
 *         con el valor de cada parámetro en cada canal de la tarjeta.
 *         Numéricos en float (unidades según la tarjeta: V, uA, V/s, s);
 *         ON/OFF y Status como entero sin signo (Status en hexadecimal).
 *
 * Compilar: gcc -Wall -DUNIX hv_readall.c -o hv_readall -lcaenhvwrapper
 *
 * Nota: verificar con
 *   grep -n "CAENHV_GetChParamInfo\|CAENHV_GetChParamProp" /usr/include/CAENHVWrapper.h
 * que las firmas coinciden con las usadas aquí (Wrapper 6.6).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "CAENHVWrapper.h"

#define MAX_CH 32

/* Cierra la sesión e informa el error de la Wrapper */
static int fail(int handle, const char *where)
{
    fprintf(stderr, "%s fallo: %s\n", where, CAENHV_GetError(handle));
    CAENHV_DeinitSystem(handle);
    return 1;
}

int main(int argc, char **argv)
{
    if (argc < 3) {
        fprintf(stderr, "Uso: %s <IP> <slot> [usuario] [password]\n", argv[0]);
        return 1;
    }
    const char *ip   = argv[1];
    unsigned short slot = (unsigned short)atoi(argv[2]);
    const char *user = (argc > 3) ? argv[3] : "admin";
    const char *pass = (argc > 4) ? argv[4] : "admin";
    int handle = -1;

    if (CAENHV_InitSystem(SY5527, LINKTYPE_TCPIP, (void *)ip,
                          user, pass, &handle) != CAENHV_OK) {
        fprintf(stderr, "InitSystem fallo: %s\n", CAENHV_GetError(handle));
        return 1;
    }

    /* Número de canales del slot desde el mapa del crate */
    unsigned short nSlots, *nChList, *serial;
    unsigned char *fwMin, *fwMax;
    char *models, *descr;
    if (CAENHV_GetCrateMap(handle, &nSlots, &nChList, &models, &descr,
                           &serial, &fwMin, &fwMax) != CAENHV_OK)
        return fail(handle, "GetCrateMap");
    if (slot >= nSlots || nChList[slot] == 0) {
        fprintf(stderr, "Slot %u vacio o inexistente\n", slot);
        CAENHV_DeinitSystem(handle);
        return 1;
    }
    unsigned short nCh = nChList[slot];
    if (nCh > MAX_CH) nCh = MAX_CH;
    CAENHV_Free(nChList); CAENHV_Free(serial); CAENHV_Free(fwMin);
    CAENHV_Free(fwMax);   CAENHV_Free(models); CAENHV_Free(descr);

    /* Lista de nombres de parámetros del canal 0 (bloques de MAX_PARAM_NAME) */
    char *parList = NULL;
    int nPar = 0;
    if (CAENHV_GetChParamInfo(handle, slot, 0, &parList, &nPar) != CAENHV_OK)
        return fail(handle, "GetChParamInfo");

    unsigned short chList[MAX_CH];
    for (unsigned short c = 0; c < nCh; c++) chList[c] = c;

    printf("Slot %u: %u canales, %d parametros\n\n", slot, nCh, nPar);

    for (int p = 0; p < nPar; p++) {
        const char *name = parList + p * MAX_PARAM_NAME;
        unsigned long type = 0;
        if (CAENHV_GetChParamProp(handle, slot, 0, name, "Type", &type) != CAENHV_OK) {
            printf("%-10s (sin tipo: %s)\n", name, CAENHV_GetError(handle));
            continue;
        }

        printf("%-10s tipo=%lu |", name, type);
        if (type == PARAM_TYPE_NUMERIC) {
            float val[MAX_CH];
            if (CAENHV_GetChParam(handle, slot, name, nCh, chList, val) != CAENHV_OK) {
                printf(" error: %s\n", CAENHV_GetError(handle));
                continue;
            }
            for (unsigned short c = 0; c < nCh; c++) printf(" %9.3f", val[c]);
        } else {
            /* ON/OFF, CHSTATUS y otros tipos enteros */
            unsigned int val[MAX_CH];
            if (CAENHV_GetChParam(handle, slot, name, nCh, chList, val) != CAENHV_OK) {
                printf(" error: %s\n", CAENHV_GetError(handle));
                continue;
            }
            for (unsigned short c = 0; c < nCh; c++)
                printf(type == PARAM_TYPE_CHSTATUS ? "    0x%04X" : " %9u", val[c]);
        }
        printf("\n");
    }

    CAENHV_Free(parList);
    CAENHV_DeinitSystem(handle);
    return 0;
}
