#!/usr/bin/env bash
# diag_sy5527.sh
# Diagnóstico de conexión entre un PC Linux (Mint/Ubuntu) y un CAEN SY5527.
# Uso:   bash diag_sy5527.sh [IP_mainframe]        (por defecto 192.168.1.117)
#        sudo bash diag_sy5527.sh ...              (para ver también el firewall)
# Salida: todo se muestra en pantalla y se guarda en diag_sy5527_<fecha>.txt
# No requiere nmap: usa ping, ip, curl y /dev/tcp de bash.

IP="${1:-192.168.1.117}"
OUT="diag_sy5527_$(date +%Y%m%d_%H%M%S).txt"
exec > >(tee "$OUT") 2>&1

sec() { echo; echo "=================== $* ==================="; }

# ---------------------------------------------------------------
sec "1. Interfaces y direcciones IPv4 de este PC"
ip -br addr
LOCAL_HIT=$(ip -o -4 addr | awk '{print $2, $4}' | grep -F " ${IP}/")
if [ -n "$LOCAL_HIT" ]; then
    echo
    echo "ALERTA: ${IP} está asignada a ESTE PC -> ${LOCAL_HIT}"
    echo "        El ping se está respondiendo a sí mismo; el mainframe no se está alcanzando."
fi

# ---------------------------------------------------------------
sec "2. Ruta hacia ${IP}"
ROUTE=$(ip route get "$IP" 2>&1)
echo "$ROUTE"
if echo "$ROUTE" | grep -q "^local"; then
    echo "ALERTA: el kernel considera ${IP} una dirección local."
fi
IFACE=$(echo "$ROUTE" | awk '{for(i=1;i<=NF;i++) if($i=="dev") print $(i+1)}' | head -n1)
echo "Interfaz de salida: ${IFACE:-desconocida}"

sec "2b. Configuración NetworkManager de las conexiones activas"
nmcli -t -f NAME,DEVICE,TYPE con show --active 2>/dev/null
for dev in $(nmcli -t -f DEVICE,TYPE dev 2>/dev/null | awk -F: '$2=="ethernet"{print $1}'); do
    echo "--- $dev"
    nmcli dev show "$dev" 2>/dev/null | grep -E "GENERAL.STATE|GENERAL.HWADDR|IP4.ADDRESS|IP4.GATEWAY"
done
echo "--- Estado físico del enlace Ethernet (carrier=1 significa cable conectado y con enlace)"
for d in /sys/class/net/*; do
    n=$(basename "$d"); [ "$n" = "lo" ] && continue
    echo "$n carrier=$(cat "$d/carrier" 2>/dev/null) operstate=$(cat "$d/operstate" 2>/dev/null)"
done

# ---------------------------------------------------------------
sec "3. Ping a ${IP}"
ping -c 3 -W 1 "$IP"
echo "(RTT < 0.1 ms suele indicar que responde el propio PC; un equipo remoto por cable suele dar > 0.1 ms)"

sec "4. Entrada ARP de ${IP}"
N=$(ip neigh show "$IP")
if [ -z "$N" ]; then
    echo "Sin entrada ARP para ${IP} (normal si la IP es local; si no, nadie respondió por Ethernet)"
else
    echo "$N"
fi

# ---------------------------------------------------------------
SUBNET="${IP%.*}"
sec "5. Barrido de ${SUBNET}.1-254 para encontrar el mainframe"
# Se guardan los PID de los pings y se espera solo a ellos.
# Un 'wait' sin argumentos esperaría también al 'tee' de la redirección y se colgaría.
PIDS=()
for i in $(seq 1 254); do
    ping -c 1 -W 1 "${SUBNET}.${i}" > /dev/null 2>&1 &
    PIDS+=($!)
done
wait "${PIDS[@]}"
echo "Equipos que respondieron por ARP (IP, interfaz, MAC):"
ip neigh show | grep -v -E "FAILED|INCOMPLETE" | grep -F "${SUBNET}." | awk '{print $1, $3, $5}'
echo "(Compara las MAC con la de la etiqueta / menú de red del SY5527.)"
echo "(Si el mainframe está en otra subred, no aparecerá aquí: leer su IP en el panel frontal.)"

# ---------------------------------------------------------------
probe_ports() {
    local host="$1"
    # 23 telnet, 80/443 web HiVoCS; 1527/4527 se prueban solo como candidatos (no verificado)
    for p in 21 22 23 80 443 1527 4527 8080; do
        if timeout 2 bash -c "echo > /dev/tcp/${host}/${p}" 2>/dev/null; then
            echo "  ${host}:${p} ABIERTO"
        else
            echo "  ${host}:${p} cerrado/sin respuesta"
        fi
    done
    echo "  HTTP:"
    curl -s -m 3 -I "http://${host}" | head -n 5 | sed 's/^/    /'
}

sec "6. Puertos TCP en ${IP}"
probe_ports "$IP"

sec "6b. Puertos en los otros equipos encontrados (excepto este PC)"
MYIPS=$(ip -o -4 addr | awk '{print $4}' | cut -d/ -f1)
for h in $(ip neigh show | grep -v -E "FAILED|INCOMPLETE" | grep -F "${SUBNET}." | awk '{print $1}'); do
    echo "$MYIPS" | grep -qx "$h" && continue
    echo "--- $h"
    probe_ports "$h"
done

# ---------------------------------------------------------------
sec "7. Firewall local"
if [ "$(id -u)" -eq 0 ]; then
    ufw status verbose 2>/dev/null || echo "ufw no instalado"
    iptables -L -n 2>/dev/null | head -n 30
else
    echo "(ejecutar con sudo para ver ufw/iptables)"
fi

# ---------------------------------------------------------------
sec "8. CAEN HV Wrapper instalada"
ldconfig -p | grep -i caenhv || echo "libcaenhvwrapper NO aparece en ldconfig"
HDR=$(find /usr /opt "$HOME" -name CAENHVWrapper.h 2>/dev/null | head -n1)
echo "Header: ${HDR:-no encontrado}"

# ---------------------------------------------------------------
sec "9. Test de conexión con la Wrapper (test_sy5527.c)"
if [ -f test_sy5527.c ] && [ -n "$HDR" ] && ldconfig -p | grep -qi caenhv; then
    gcc -Wall -I"$(dirname "$HDR")" test_sy5527.c -o test_sy5527 -lcaenhvwrapper \
        && ./test_sy5527 "$IP" admin admin
else
    echo "Se omite: falta test_sy5527.c en este directorio, el header o la librería."
fi

echo
echo "Reporte guardado en: $OUT"
