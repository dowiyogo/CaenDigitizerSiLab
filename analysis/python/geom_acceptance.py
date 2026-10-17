#!/usr/bin/env python3
"""
geom_acceptance.py
Monte Carlo geométrico de aceptancia del trigger vs hodoscopios.
NO requiere Geant4 — solo trazado de rayos rectos con distribución cos²θ.

Geometría (z creciente hacia arriba, muones descienden):
  hod2018  z=0 mm            área activa ±16.5 mm
  hod2019  z=98.4 mm         área activa ±16.5 mm
  T1 center z=124.1 mm       150×150 mm, centro en (-42.5, -42.5) mm
  T2 center z=134.1 mm       150×150 mm, centro en (+42.5, +42.5) mm
  Muones generados z=160 mm y descienden.
"""

import numpy as np
import os
import sys

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False
    print("WARNING: matplotlib not available; figures will not be generated")

# ─── Geometría ───────────────────────────────────────────────────────────────
Z_START = 160.0   # mm — generación de muones
Z_T2    = 134.1   # mm — T2 (primera placa que ven los muones al descender)
Z_T1    = 124.1   # mm — T1
Z_HOD19 =  98.4   # mm — hod2019
Z_HOD18 =   0.0   # mm — hod2018

# T1: centro (-42.5, -42.5), 150×150
T1_xlo, T1_xhi = -42.5-75, -42.5+75  # [-117.5,  32.5]
T1_ylo, T1_yhi = -42.5-75, -42.5+75

# T2: centro (+42.5, +42.5), 150×150
T2_xlo, T2_xhi =  42.5-75,  42.5+75  # [ -32.5, 117.5]
T2_ylo, T2_yhi =  42.5-75,  42.5+75

HOD_HALF = 16.5   # ±16.5 mm → 33×33 mm área activa

# ─── Generación ──────────────────────────────────────────────────────────────
def run_mc(N=1_000_000, seed=42, outdir="."):
    rng = np.random.default_rng(seed)
    print(f"Generando {N:,} muones (seed={seed})...")

    # Ángulo cenital: distribución cos²θ sinθ
    # CDF: F(θ)=1−cos³θ → θ = arccos(U^{1/3})
    U   = rng.uniform(0, 1, N)
    theta = np.arccos(U**(1.0/3.0))
    phi   = rng.uniform(0, 2*np.pi, N)

    # Vectores de dirección (descenso en -z)
    dx = np.sin(theta) * np.cos(phi)
    dy = np.sin(theta) * np.sin(phi)
    dz = -np.cos(theta)   # siempre negativo

    # Posición inicial (amplia para capturar todos los rayos que pueden cruzar el trigger)
    # El overlap del trigger es 65×65 centrado en (0,0) a z≈129. Con θ_max≈45° y
    # dz_vertical≈31 mm, el origen puede estar hasta ±75+31≈106 mm del centro → usar ±180 mm.
    x0 = rng.uniform(-180, 180, N)
    y0 = rng.uniform(-180, 180, N)
    z0 = Z_START

    # Propagación a cada plano: t = (z_target - z0) / dz
    def prop(z_target):
        t = (z_target - z0) / dz
        return x0 + dx*t, y0 + dy*t

    xT2, yT2   = prop(Z_T2)
    xT1, yT1   = prop(Z_T1)
    x19, y19   = prop(Z_HOD19)
    x18, y18   = prop(Z_HOD18)

    # Condiciones de trigger (geométricas)
    inT2 = (xT2 >= T2_xlo) & (xT2 <= T2_xhi) & (yT2 >= T2_ylo) & (yT2 <= T2_yhi)
    inT1 = (xT1 >= T1_xlo) & (xT1 <= T1_xhi) & (yT1 >= T1_ylo) & (yT1 <= T1_yhi)
    trig = inT1 & inT2

    # Áreas activas de hodoscopios
    in18 = (np.abs(x18) <= HOD_HALF) & (np.abs(y18) <= HOD_HALF)
    in19 = (np.abs(x19) <= HOD_HALF) & (np.abs(y19) <= HOD_HALF)

    N_trig = trig.sum()
    if N_trig == 0:
        print("ERROR: no triggered muons!")
        return

    f18   = (trig & in18).sum() / N_trig
    f19   = (trig & in19).sum() / N_trig
    fboth = (trig & in18 & in19).sum() / N_trig

    # Ángulos para los eventos triggerados (para fase-space plot)
    theta_trig = theta[trig]
    theta_both = theta[trig & in18 & in19]

    # ── Resultados ────────────────────────────────────────────────────────────
    err18   = np.sqrt(f18   * (1-f18)   / N_trig)
    err19   = np.sqrt(f19   * (1-f19)   / N_trig)
    errboth = np.sqrt(fboth * (1-fboth) / N_trig)

    print(f"\n=== RESULTADOS DE ACEPTANCIA GEOMÉTRICA ===")
    print(f"  N generados:   {N:,}")
    print(f"  N trigger:     {N_trig:,}  ({100*N_trig/N:.3f}% del total)")
    print(f"  f_2018         = {f18:.4f} ± {err18:.4f}   ({100*f18:.2f}%)")
    print(f"  f_2019         = {f19:.4f} ± {err19:.4f}   ({100*f19:.2f}%)")
    print(f"  f_both         = {fboth:.4f} ± {errboth:.4f}   ({100*fboth:.2f}%)")
    print(f"\n  Implicación: solo {100*fboth:.1f}% de los muones triggerados cruzan ambas")
    print(f"  áreas activas. ε_det(hod2018) = ε_sig / {f18:.3f}")

    # ── Reconciliación de eficiencias ─────────────────────────────────────────
    # Valores del análisis experimental
    eps_sys_v7  = 0.999    # fracción con cualquier hit en hod2018 (v7)
    eps_sig_v7  = 0.1227   # fracción con cluster limpio Q>umbral (v7, charge_analysis)
    eps_sig_v8  = 0.0249   # v8
    eps_det_v7  = eps_sig_v7  / f18 if f18 > 0 else 0
    eps_det_v8  = eps_sig_v8  / f18 if f18 > 0 else 0

    print(f"\n=== RECONCILIACIÓN DE EFICIENCIAS hod2018X ===")
    print(f"  ε_sys   (v7) = {eps_sys_v7*100:.1f}%  [fondo incluido — hits de cualquier canal]")
    print(f"  ε_sig   (v7) = {eps_sig_v7*100:.2f}%  [cluster limpio Q>umbral_Q, cualquier posición]")
    print(f"  ε_det   (v7) = {eps_det_v7*100:.2f}%  [ε_sig / f_2018 — eficiencia de detector corregida]")
    print(f"  ε_sig   (v8) = {eps_sig_v8*100:.2f}%")
    print(f"  ε_det   (v8) = {eps_det_v8*100:.2f}%")
    print(f"\n  CAVEAT: ε_sys=99.9% NO es eficiencia de muones.")
    print(f"  Es tasa de sistema (include fondo). ε_det es la estimación corregida.")

    # ── Guardar resultados ────────────────────────────────────────────────────
    with open(os.path.join(outdir, "geom_acceptance.txt"), "w") as fout:
        fout.write("=== Aceptancia Geométrica — Monte Carlo ===\n\n")
        fout.write(f"N_generados   = {N:,}\n")
        fout.write(f"N_trigger     = {N_trig:,}  ({100*N_trig/N:.3f}%)\n")
        fout.write(f"f_2018        = {f18:.5f} ± {err18:.5f}  ({100*f18:.2f}%)\n")
        fout.write(f"f_2019        = {f19:.5f} ± {err19:.5f}  ({100*f19:.2f}%)\n")
        fout.write(f"f_both        = {fboth:.5f} ± {errboth:.5f}  ({100*fboth:.2f}%)\n\n")
        fout.write("Reconciliación eficiencias hod2018X:\n")
        fout.write(f"  eps_sys_v7  = {eps_sys_v7*100:.1f}%  [tasa de sistema, fondo incluido]\n")
        fout.write(f"  eps_sig_v7  = {eps_sig_v7*100:.2f}%  [cluster limpio Q>umbral]\n")
        fout.write(f"  eps_det_v7  = {eps_det_v7*100:.2f}%  [= eps_sig / f_2018, corregida por geometria]\n")
        fout.write(f"  eps_sig_v8  = {eps_sig_v8*100:.2f}%\n")
        fout.write(f"  eps_det_v8  = {eps_det_v8*100:.2f}%\n")
    print(f"\n  → geom_acceptance.txt escrito en {outdir}")

    if not HAS_MPL:
        print("Skipping figures (matplotlib unavailable)")
        return f18, f19, fboth

    # ── Figura 1: Vista superior — impactos en hod2018 de muones triggerados ──
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle("Aceptancia geométrica: trigger vs hodoscopios", fontsize=13, fontweight='bold')

    # Panel izquierdo: posiciones en hod2018 para muones trigger
    ax = axes[0]
    sample = np.random.choice(N_trig, min(N_trig, 30000), replace=False)
    ax.scatter(x18[trig][sample], y18[trig][sample], s=0.3, alpha=0.3, color='steelblue',
               label=f'Trigger ({100*N_trig/N:.2f}% del total)')
    # Área activa
    rect = plt.Rectangle((-HOD_HALF, -HOD_HALF), 2*HOD_HALF, 2*HOD_HALF,
                          fill=False, edgecolor='red', linewidth=2,
                          label=f'Área activa hod2018 (33×33 mm)')
    ax.add_patch(rect)
    # Overlap 65×65
    rect2 = plt.Rectangle((-32.5, -32.5), 65, 65,
                           fill=False, edgecolor='orange', linewidth=1.5,
                           linestyle='--', label='Overlap T1∧T2 (65×65 mm)')
    ax.add_patch(rect2)
    ax.set_xlim(-100, 100); ax.set_ylim(-100, 100)
    ax.set_xlabel("X en hod2018 [mm]"); ax.set_ylabel("Y en hod2018 [mm]")
    ax.set_title(f"Vista superior: impactos en hod2018\n"
                 f"f_2018={100*f18:.1f}%  f_both={100*fboth:.1f}%", fontsize=11)
    ax.legend(fontsize=8, loc='upper right'); ax.set_aspect('equal')
    ax.text(0.02, 0.02, f"N_trigger={N_trig:,}", transform=ax.transAxes, fontsize=9)

    # Panel derecho: distribución de θ — trigger vs ambos hodoscopios
    ax = axes[1]
    bins = np.linspace(0, 45, 50)
    ax.hist(np.degrees(theta_trig), bins=bins, density=True, alpha=0.6, color='steelblue',
            label=f'Trigger ({N_trig:,})')
    ax.hist(np.degrees(theta_both), bins=bins, density=True, alpha=0.6, color='firebrick',
            label=f'Trigger ∧ hod18 ∧ hod19 ({(trig & in18 & in19).sum():,})')
    # cos²θ reference
    th_plot = np.linspace(0, 45, 200)
    cos2_dist = 3 * np.cos(np.radians(th_plot))**2 * np.sin(np.radians(th_plot)) * np.pi/180
    ax.plot(th_plot, cos2_dist / cos2_dist.max() * ax.get_ylim()[1] if ax.get_ylim()[1]>0 else cos2_dist,
            'k--', lw=1.5, label='cos²θ sinθ (teórico)', alpha=0.7)
    ax.set_xlabel("Ángulo cenital θ [grados]"); ax.set_ylabel("Densidad")
    ax.set_title("Fase espacial angular\n(muones de mayor θ → menos probabilidad de cruzar hod activo)",
                 fontsize=10)
    ax.legend(fontsize=9)

    plt.tight_layout()
    fig.savefig(os.path.join(outdir, "acceptance_MC_topview.png"), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("  → acceptance_MC_topview.png generada")

    # ── Figura 2: Fase espacial (2D) ─────────────────────────────────────────
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle("Distribuciones angulares: trigger vs telescopio completo", fontsize=12)

    # θ vs φ para trigger
    ax = axes[0]
    h, xe, ye = np.histogram2d(np.degrees(theta[trig]), np.degrees(phi[trig]), bins=40)
    pcm = ax.pcolormesh(ye, xe, h, cmap='Blues')
    fig.colorbar(pcm, ax=ax, label='Counts')
    ax.set_xlabel("φ [grados]"); ax.set_ylabel("θ [grados]")
    ax.set_title(f"Trigger (T1∧T2)\nN={N_trig:,}")

    # θ vs φ para ambos hodoscopios
    ax = axes[1]
    both_mask = trig & in18 & in19
    h2, xe, ye = np.histogram2d(np.degrees(theta[both_mask]), np.degrees(phi[both_mask]), bins=40)
    pcm2 = ax.pcolormesh(ye, xe, h2, cmap='Reds')
    fig.colorbar(pcm2, ax=ax, label='Counts')
    ax.set_xlabel("φ [grados]"); ax.set_ylabel("θ [grados]")
    ax.set_title(f"Trigger ∧ hod2018 ∧ hod2019\nN={(both_mask).sum():,}  (f_both={100*fboth:.1f}%)")

    plt.tight_layout()
    fig.savefig(os.path.join(outdir, "acceptance_phase_space.png"), dpi=150, bbox_inches='tight')
    plt.close(fig)
    print("  → acceptance_phase_space.png generada")

    return f18, f19, fboth

# ─── Main ────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    import datetime
    ds = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    outdir = os.path.join("outputs", f"geom_acceptance_{ds}")
    os.makedirs(outdir, exist_ok=True)
    print(f"=== geom_acceptance.py ===")
    print(f"Outdir: {outdir}\n")
    f18, f19, fboth = run_mc(N=1_000_000, seed=42, outdir=outdir)
    print(f"\nDone. Output: {outdir}")
