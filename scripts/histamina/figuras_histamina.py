"""Figuras 35 y 36 y estadísticos de histamina en recepción.

Lee la hoja Base_Normalizada de Jun_Control_Calidad_Histamina_v2.xlsx
(sin dependencias externas para el libro) y escribe:
  figuras/histamina/fig35_histamina_especie.png
  figuras/histamina/fig36_histamina_imr.png
  figuras/histamina/estadisticos_histamina.json

Uso:
  python3 scripts/histamina/figuras_histamina.py [ruta_al_xlsx]
Requiere matplotlib y numpy.
"""
import json
import math
import re
import statistics as st
import sys
import zipfile
import xml.etree.ElementTree as ET
from collections import Counter, defaultdict
from datetime import date, timedelta
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RAIZ = Path(__file__).resolve().parents[2]
XLSX = Path(sys.argv[1]) if len(sys.argv) > 1 else RAIZ / "Jun_Control_Calidad_Histamina_v2.xlsx"
SALIDA = RAIZ / "figuras" / "histamina"
NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"

# Constantes de cartas I-MR para n = 2 (AIAG/Montgomery)
D4 = 3.267
D2 = 1.128


def leer_hoja(xlsx, nombre_archivo_hoja):
    """Devuelve las filas de una hoja como listas de texto, sin cabecera de libro."""
    z = zipfile.ZipFile(xlsx)
    raiz = ET.fromstring(z.read(nombre_archivo_hoja))
    filas = []
    for fila in raiz.iter(NS + "row"):
        celdas = {}
        for c in fila.findall(NS + "c"):
            ref = re.match(r"[A-Z]+", c.get("r")).group(0)
            n = 0
            for ch in ref:
                n = n * 26 + ord(ch) - 64
            t = c.get("t")
            v = c.find(NS + "v")
            if t == "inlineStr":
                val = "".join(x.text or "" for x in c.iter(NS + "t"))
            else:
                val = "" if v is None else v.text
            celdas[n - 1] = val
        if celdas:
            filas.append([celdas.get(k, "") for k in range(max(celdas) + 1)])
    return filas


def cargar_registros(xlsx):
    filas = leer_hoja(xlsx, "xl/worksheets/sheet4.xml")  # Base_Normalizada
    cab = next(i for i, f in enumerate(filas) if f and f[0] == "ID")
    registros = []
    for f in filas[cab + 1:]:
        if len(f) < 12 or not f[0].strip():
            continue
        registros.append({
            "id": int(float(f[0])),
            "fecha": date(1899, 12, 30) + timedelta(days=int(float(f[1]))),
            "lote": f[5],
            "especie": f[6],
            "estado": f[11],
            "ppm": float(f[9]),
        })
    return registros


def mann_whitney_z(a, b):
    """Z de Mann-Whitney con corrección por empates; devuelve (z, p bilateral)."""
    datos = sorted([(x, 0) for x in a] + [(x, 1) for x in b])
    rangos, i, tie = {}, 0, 0
    while i < len(datos):
        j = i
        while j < len(datos) and datos[j][0] == datos[i][0]:
            j += 1
        r = (i + j + 1) / 2
        for k in range(i, j):
            rangos[k] = r
        t = j - i
        tie += t ** 3 - t
        i = j
    r2 = sum(rangos[k] for k in range(len(datos)) if datos[k][1] == 1)
    n1, n2 = len(a), len(b)
    u2 = r2 - n2 * (n2 + 1) / 2
    u = min(n1 * n2 - u2, u2)
    big_n = n1 + n2
    sigma = math.sqrt(n1 * n2 / 12 * ((big_n + 1) - tie / (big_n * (big_n - 1))))
    z = (u - n1 * n2 / 2) / sigma
    return z, math.erfc(abs(z) / math.sqrt(2))


def kruskal_wallis(grupos):
    """H con corrección por empates, eta² de Kruskal-Wallis (Tomczak y Tomczak)."""
    todos = sorted((v, g) for g, vals in grupos.items() for v in vals)
    n = len(todos)
    suma = defaultdict(float)
    i, tie = 0, 0
    while i < n:
        j = i
        while j < n and todos[j][0] == todos[i][0]:
            j += 1
        r = (i + j + 1) / 2
        t = j - i
        tie += t ** 3 - t
        for k in range(i, j):
            suma[todos[k][1]] += r
        i = j
    h = 12 / (n * (n + 1)) * sum(suma[g] ** 2 / len(v) for g, v in grupos.items()) - 3 * (n + 1)
    h = h / (1 - tie / (n ** 3 - n))
    k = len(grupos)
    eta2 = (h - k + 1) / (n - k)
    return h, k - 1, eta2


def estadisticos(regs):
    valores = [r["ppm"] for r in regs]
    por_especie = defaultdict(list)
    for r in regs:
        por_especie[r["especie"]].append(r["ppm"])
    alb = [r["ppm"] for r in regs if r["especie"] == "ALBACORA"]
    bar_todas = [r["ppm"] for r in regs if r["especie"] == "SKJ"]
    bar_validas = [r["ppm"] for r in regs if r["especie"] == "SKJ" and r["estado"] == "Válido"]
    z_todas, p_todas = mann_whitney_z(alb, bar_todas)
    z_val, p_val = mann_whitney_z(alb, bar_validas)
    lotes = defaultdict(list)
    for r in regs:
        lotes[r["lote"]].append(r["ppm"])
    h, gl, eta2 = kruskal_wallis(lotes)

    # Carta I-MR sobre lecturas en orden de identificador
    orden = sorted(regs, key=lambda r: r["id"])
    x = np.array([r["ppm"] for r in orden])
    mr = np.abs(np.diff(x))
    xbar, mrbar = x.mean(), mr.mean()
    sigma = mrbar / D2
    lsc, lic = xbar + 3 * sigma, xbar - 3 * sigma
    lsc_mr = D4 * mrbar
    fuera_i = int(((x > lsc) | (x < lic)).sum())
    fuera_mr = int((mr > lsc_mr).sum())

    def resumen(v):
        v = list(v)
        return {
            "n": len(v),
            "mediana": st.median(v),
            "media": st.mean(v),
            "de": st.stdev(v),
            "min": min(v),
            "max": max(v),
            "pct_ge_10": 100 * sum(1 for a in v if a >= 10) / len(v),
        }

    return {
        "global": resumen(valores),
        "pct_ge_5": 100 * sum(1 for a in valores if a >= 5) / len(valores),
        "n_ge_5": sum(1 for a in valores if a >= 5),
        "especies": {e: resumen(v) for e, v in por_especie.items()},
        "mwu_albacora_vs_skj_todas": {"z": z_todas, "p": p_todas},
        "mwu_albacora_vs_skj_validas": {"z": z_val, "p": p_val},
        "kruskal_wallis_lotes": {"H": h, "gl": gl, "eta2_H": eta2, "n_lotes": len(lotes)},
        "imr": {
            "xbar": xbar, "mrbar": mrbar, "sigma": sigma,
            "lsc_i": lsc, "lic_i": lic, "lsc_mr": lsc_mr,
            "fuera_i": fuera_i, "fuera_mr": fuera_mr, "n": len(x),
        },
    }, orden


def figura_35(regs, est, salida):
    grupos = [
        ("ALBACORA", "Albacora"),
        ("SKJ", "Barrilete"),
        ("SIN ESPECIFICAR", "Sin especie declarada"),
    ]
    fig, ax = plt.subplots(figsize=(7.2, 4.8), dpi=200)
    datos = []
    for i, (clave, _) in enumerate(grupos, start=1):
        datos.append([r["ppm"] for r in regs if r["especie"] == clave])
    ax.boxplot(datos, positions=range(1, 4), widths=0.45, showfliers=False,
               medianprops=dict(color="black", linewidth=1.4),
               boxprops=dict(color="#555555"), whiskerprops=dict(color="#555555"),
               capprops=dict(color="#555555"))
    rng = np.random.default_rng(35)
    for i, (clave, _) in enumerate(grupos, start=1):
        sel = [r for r in regs if r["especie"] == clave]
        xs = i + rng.uniform(-0.18, 0.18, size=len(sel))
        for x_, r in zip(xs, sel):
            if clave == "SKJ" and r["estado"] != "Válido":
                ax.scatter(x_, r["ppm"], s=14, facecolors="none", edgecolors="#B5651D", linewidths=0.8, zorder=3)
            else:
                color = {"ALBACORA": "#1F4E79", "SKJ": "#B5651D", "SIN ESPECIFICAR": "#7F7F7F"}[clave]
                ax.scatter(x_, r["ppm"], s=14, color=color, alpha=0.75, zorder=3)
    ax.axhline(10, color="#C00000", linestyle="--", linewidth=1)
    ax.axhline(30, color="#7030A0", linestyle=":", linewidth=1.2)
    ax.text(0.55, 10.4, "10 ppm: umbral de tamizaje", color="#C00000", fontsize=7.5, ha="left", va="bottom")
    ax.text(0.55, 30.4, "30 ppm: límite por pez individual", color="#7030A0", fontsize=7.5, ha="left", va="bottom")
    for i, (clave, _) in enumerate(grupos, start=1):
        e = est["especies"][clave]
        ax.text(i, -4.6, f"n = {e['n']}\nmed. {e['mediana']:.2f}\nmedia {e['media']:.2f}",
                ha="center", va="top", fontsize=7.5)
    ax.set_xticks(range(1, 4))
    ax.set_xticklabels([etq for _, etq in grupos], fontsize=9)
    ax.set_ylabel("Histamina (ppm)")
    ax.set_ylim(-9, 38)
    ax.set_yticks(range(0, 35, 5))
    ax.grid(axis="y", color="#DDDDDD", linewidth=0.6)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    m1 = est["mwu_albacora_vs_skj_todas"]
    m2 = est["mwu_albacora_vs_skj_validas"]
    kw = est["kruskal_wallis_lotes"]
    texto = (f"Mann-Whitney albacora vs. barrilete: z = {m1['z']:.2f} (p < 0.001); "
             f"sin duplicados: z = {m2['z']:.2f} (p < 0.0001)\n"
             f"Kruskal-Wallis entre {kw['n_lotes']} lotes: H = {kw['H']:.1f}, gl = {kw['gl']}, "
             f"p < 0.0001, η²_H = {kw['eta2_H']:.3f}")
    ax.text(0.01, 0.985, texto, transform=ax.transAxes, fontsize=6.8, va="top",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#BBBBBB"))
    fig.tight_layout()
    fig.savefig(salida / "fig35_histamina_especie.png")
    plt.close(fig)


def figura_36(orden, est, salida):
    x = np.array([r["ppm"] for r in orden])
    idx = np.arange(1, len(x) + 1)
    mr = np.abs(np.diff(x))
    im = est["imr"]
    fig, (a1, a2) = plt.subplots(2, 1, figsize=(7.2, 5.6), dpi=200, sharex=True,
                                 gridspec_kw={"height_ratios": [2.2, 1]})
    a1.plot(idx, x, color="#1F4E79", linewidth=0.7, marker="o", markersize=2.2)
    a1.axhline(im["xbar"], color="black", linewidth=1, label=f"LC = {im['xbar']:.2f} ppm")
    a1.axhline(im["lsc_i"], color="#C00000", linestyle="--", linewidth=1,
               label=f"LSC = {im['lsc_i']:.2f} ppm")
    a1.axhline(max(im["lic_i"], 0), color="#C00000", linestyle="--", linewidth=1,
               label=f"LIC = {max(im['lic_i'], 0):.2f} ppm (truncado en 0)")
    a1.set_ylabel("Histamina (ppm)")
    a1.legend(fontsize=7, loc="upper right", frameon=False)
    a1.grid(axis="y", color="#DDDDDD", linewidth=0.6)
    a2.plot(idx[1:], mr, color="#555555", linewidth=0.7, marker="o", markersize=2)
    a2.axhline(im["mrbar"], color="black", linewidth=1)
    a2.axhline(im["lsc_mr"], color="#C00000", linestyle="--", linewidth=1)
    a2.set_ylabel("MR (ppm)")
    a2.set_xlabel("Lectura de rutina (orden de identificador, 1 a 524)")
    a2.grid(axis="y", color="#DDDDDD", linewidth=0.6)
    for a in (a1, a2):
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
    texto = (f"n = {im['n']}; MR media = {im['mrbar']:.3f} ppm; sigma estimada = MR media / d2 = {im['sigma']:.3f} ppm; "
             f"lecturas fuera de ±3 sigma: {im['fuera_i']}; MR fuera de LSC: {im['fuera_mr']}")
    fig.text(0.01, 0.005, texto, fontsize=6.8)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(salida / "fig36_histamina_imr.png")
    plt.close(fig)


def main():
    SALIDA.mkdir(parents=True, exist_ok=True)
    regs = cargar_registros(XLSX)
    assert len(regs) == 524, f"se esperaban 524 registros, hay {len(regs)}"
    est, orden = estadisticos(regs)
    figura_35(regs, est, SALIDA)
    figura_36(orden, est, SALIDA)
    (SALIDA / "estadisticos_histamina.json").write_text(
        json.dumps(est, indent=2, ensure_ascii=False, default=float), encoding="utf8")
    print(json.dumps(est, indent=2, ensure_ascii=False, default=float))


if __name__ == "__main__":
    main()
