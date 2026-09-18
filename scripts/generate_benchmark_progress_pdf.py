#!/usr/bin/env python3
"""PDF de avance de CHA-MD-BA, con énfasis en el benchmark 6M03."""

from __future__ import annotations

import html
from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import cm, mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    HRFlowable,
    KeepTogether,
    ListFlowable,
    ListItem,
    PageBreak,
    Paragraph,
    Preformatted,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)

ROOT = Path(__file__).resolve().parents[1]
FONT_DIR = Path("/usr/share/fonts/truetype/dejavu")
OUTPUT = ROOT / "docs" / "benchmarks" / "CHA-MD-BA_Benchmark_6M03_avance.pdf"

NAVY = colors.HexColor("#1B365D")
TEAL = colors.HexColor("#1F6F7A")
GOLD = colors.HexColor("#C4A35A")
INK = colors.HexColor("#1A1A1A")
MUTED = colors.HexColor("#4A5560")
RULE = colors.HexColor("#D6DDE5")
ROW = colors.HexColor("#F4F7FA")
HEAD_BG = colors.HexColor("#1B365D")
CODE_BG = colors.HexColor("#F3F1EB")
WARN_BG = colors.HexColor("#F8F1E3")


def register_fonts() -> None:
    pdfmetrics.registerFont(TTFont("DejaVu", str(FONT_DIR / "DejaVuSans.ttf")))
    pdfmetrics.registerFont(TTFont("DejaVu-Bold", str(FONT_DIR / "DejaVuSans-Bold.ttf")))
    pdfmetrics.registerFont(TTFont("DejaVu-Italic", str(FONT_DIR / "DejaVuSans.ttf")))
    pdfmetrics.registerFont(TTFont("DejaVuMono", str(FONT_DIR / "DejaVuSansMono.ttf")))
    pdfmetrics.registerFontFamily(
        "DejaVu",
        normal="DejaVu",
        bold="DejaVu-Bold",
        italic="DejaVu",
        boldItalic="DejaVu-Bold",
    )


def styles():
    s = {
        "kicker": ParagraphStyle(
            "kicker",
            fontName="DejaVu-Bold",
            fontSize=9,
            textColor=TEAL,
            letterSpacing=1.2,
            alignment=TA_CENTER,
            spaceAfter=6,
        ),
        "title": ParagraphStyle(
            "title",
            fontName="DejaVu-Bold",
            fontSize=22,
            leading=26,
            textColor=NAVY,
            alignment=TA_CENTER,
            spaceAfter=8,
        ),
        "subtitle": ParagraphStyle(
            "subtitle",
            fontName="DejaVu",
            fontSize=12,
            leading=16,
            textColor=MUTED,
            alignment=TA_CENTER,
            spaceAfter=4,
        ),
        "meta": ParagraphStyle(
            "meta",
            fontName="DejaVu",
            fontSize=9.5,
            leading=13,
            textColor=MUTED,
            alignment=TA_CENTER,
            spaceAfter=2,
        ),
        "h1": ParagraphStyle(
            "h1",
            fontName="DejaVu-Bold",
            fontSize=13.5,
            leading=17,
            textColor=NAVY,
            spaceBefore=16,
            spaceAfter=8,
        ),
        "h2": ParagraphStyle(
            "h2",
            fontName="DejaVu-Bold",
            fontSize=11.5,
            leading=15,
            textColor=TEAL,
            spaceBefore=12,
            spaceAfter=6,
        ),
        "body": ParagraphStyle(
            "body",
            fontName="DejaVu",
            fontSize=10,
            leading=14,
            textColor=INK,
            alignment=TA_JUSTIFY,
            spaceAfter=8,
        ),
        "bullet": ParagraphStyle(
            "bullet",
            fontName="DejaVu",
            fontSize=10,
            leading=13.5,
            textColor=INK,
            leftIndent=4,
        ),
        "cell": ParagraphStyle(
            "cell",
            fontName="DejaVu",
            fontSize=8.5,
            leading=11.5,
            textColor=INK,
        ),
        "cell_h": ParagraphStyle(
            "cell_h",
            fontName="DejaVu-Bold",
            fontSize=8.5,
            leading=11.5,
            textColor=colors.white,
        ),
        "caption": ParagraphStyle(
            "caption",
            fontName="DejaVu",
            fontSize=8,
            leading=11,
            textColor=MUTED,
            spaceBefore=3,
            spaceAfter=10,
            alignment=TA_LEFT,
        ),
        "code": ParagraphStyle(
            "code",
            fontName="DejaVuMono",
            fontSize=8,
            leading=11,
            textColor=INK,
            backColor=CODE_BG,
            leftIndent=6,
            rightIndent=6,
            spaceBefore=4,
            spaceAfter=10,
        ),
        "callout": ParagraphStyle(
            "callout",
            fontName="DejaVu",
            fontSize=9.5,
            leading=13,
            textColor=INK,
            alignment=TA_LEFT,
        ),
        "footer": ParagraphStyle(
            "footer",
            fontName="DejaVu",
            fontSize=8,
            textColor=MUTED,
        ),
        "cover_small": ParagraphStyle(
            "cover_small",
            fontName="DejaVu",
            fontSize=10,
            leading=14,
            textColor=INK,
            alignment=TA_CENTER,
        ),
    }
    return s


def md(text: str) -> str:
    """Escape HTML and turn **bold** / `code` into reportlab markup."""
    text = html.escape(text)
    while "**" in text:
        start = text.find("**")
        end = text.find("**", start + 2)
        if end == -1:
            break
        text = text[:start] + "<b>" + text[start + 2 : end] + "</b>" + text[end + 2 :]
    while "`" in text:
        start = text.find("`")
        end = text.find("`", start + 1)
        if end == -1:
            break
        text = (
            text[:start]
            + '<font name="DejaVuMono" size="8.5" color="#1F6F7A">'
            + text[start + 1 : end]
            + "</font>"
            + text[end + 1 :]
        )
    return text


def P(text: str, style) -> Paragraph:
    return Paragraph(md(text), style)


def bullets(items, st) -> ListFlowable:
    return ListFlowable(
        [ListItem(P(item, st["bullet"]), leftIndent=12, bulletColor=TEAL) for item in items],
        bulletType="bullet",
        start="•",
        leftIndent=18,
        bulletFontName="DejaVu",
        bulletFontSize=10,
        spaceBefore=2,
        spaceAfter=8,
    )


def table(headers, rows, col_widths, st) -> Table:
    data = [[Paragraph(md(h), st["cell_h"]) for h in headers]]
    for row in rows:
        data.append([Paragraph(md(c), st["cell"]) for c in row])
    t = Table(data, colWidths=col_widths, repeatRows=1)
    style_cmds = [
        ("BACKGROUND", (0, 0), (-1, 0), HEAD_BG),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), "DejaVu-Bold"),
        ("FONTNAME", (0, 1), (-1, -1), "DejaVu"),
        ("FONTSIZE", (0, 0), (-1, -1), 8.5),
        ("ALIGN", (0, 0), (-1, 0), "LEFT"),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("LEFTPADDING", (0, 0), (-1, -1), 6),
        ("RIGHTPADDING", (0, 0), (-1, -1), 6),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
        ("GRID", (0, 0), (-1, -1), 0.3, RULE),
        ("LINEBELOW", (0, 0), (-1, 0), 1.2, GOLD),
    ]
    for i in range(1, len(data)):
        if i % 2 == 0:
            style_cmds.append(("BACKGROUND", (0, i), (-1, i), ROW))
    t.setStyle(TableStyle(style_cmds))
    return t


def callout(title: str, body: str, st) -> Table:
    inner = [
        [Paragraph(md(f"**{title}**"), ParagraphStyle("ct", parent=st["callout"], fontName="DejaVu-Bold", textColor=NAVY))],
        [Paragraph(md(body), st["callout"])],
    ]
    box = Table(inner, colWidths=[16.5 * cm])
    box.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, -1), WARN_BG),
                ("BOX", (0, 0), (-1, -1), 0.6, GOLD),
                ("LEFTPADDING", (0, 0), (-1, -1), 10),
                ("RIGHTPADDING", (0, 0), (-1, -1), 10),
                ("TOPPADDING", (0, 0), (0, 0), 8),
                ("BOTTOMPADDING", (0, -1), (-1, -1), 8),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
            ]
        )
    )
    return box


def rule() -> HRFlowable:
    return HRFlowable(width="100%", thickness=1.1, color=GOLD, spaceBefore=4, spaceAfter=10)


def header_footer(canvas, doc):
    canvas.saveState()
    w, h = A4
    canvas.setFillColor(NAVY)
    canvas.rect(0, h - 12 * mm, w, 12 * mm, fill=1, stroke=0)
    canvas.setFillColor(GOLD)
    canvas.rect(0, h - 12.8 * mm, w, 1.6 * mm, fill=1, stroke=0)
    canvas.setFillColor(colors.white)
    canvas.setFont("DejaVu", 8)
    canvas.drawString(18 * mm, h - 8 * mm, "CHA-MD-BA  ·  Informe de avance")
    canvas.drawRightString(w - 18 * mm, h - 8 * mm, "Benchmark 6M03 · NaCl 0.15 M · 310 K")

    canvas.setFillColor(NAVY)
    canvas.rect(0, 0, w, 12 * mm, fill=1, stroke=0)
    canvas.setFillColor(GOLD)
    canvas.rect(0, 12 * mm, w, 1.2 * mm, fill=1, stroke=0)
    canvas.setFillColor(colors.white)
    canvas.setFont("DejaVu", 8)
    canvas.drawString(18 * mm, 5 * mm, "github.com/mixcoha/cha-MD-ba")
    canvas.drawRightString(w - 18 * mm, 5 * mm, f"pág. {doc.page}")
    canvas.restoreState()


def cover_page(canvas, doc):
    header_footer(canvas, doc)
    w, h = A4
    canvas.saveState()
    canvas.setFillColor(NAVY)
    canvas.rect(0, h - 108 * mm, w, 96 * mm, fill=1, stroke=0)
    canvas.setFillColor(GOLD)
    canvas.rect(0, h - 110 * mm, w, 2.4 * mm, fill=1, stroke=0)
    canvas.setFillColor(GOLD)
    canvas.setFont("DejaVu-Bold", 9)
    canvas.drawCentredString(w / 2.0, h - 40 * mm, "PIPELINE DE DINÁMICA MOLECULAR")
    canvas.setFillColor(colors.white)
    canvas.setFont("DejaVu-Bold", 28)
    canvas.drawCentredString(w / 2.0, h - 56 * mm, "CHA-MD-BA")
    canvas.setFont("DejaVu", 11.5)
    canvas.drawCentredString(w / 2.0, h - 72 * mm, "Informe de avance del benchmark")
    canvas.drawCentredString(w / 2.0, h - 79 * mm, "SARS-CoV-2 Mpro apo  ·  PDB 6M03")
    canvas.setFillColor(colors.HexColor("#D6DDE5"))
    canvas.setFont("DejaVu", 9.5)
    canvas.drawCentredString(
        w / 2.0, h - 94 * mm, "NaCl 0.15 M   ·   310 K   ·   1 bar   ·   10 ns de referencia"
    )
    canvas.restoreState()


def build(output: Path = OUTPUT) -> Path:
    register_fonts()
    st = styles()
    output.parent.mkdir(parents=True, exist_ok=True)

    doc = SimpleDocTemplate(
        str(output),
        pagesize=A4,
        leftMargin=18 * mm,
        rightMargin=18 * mm,
        topMargin=20 * mm,
        bottomMargin=18 * mm,
        title="CHA-MD-BA — Avance del benchmark 6M03",
        author="Edgar Mixcoha",
        subject="Informe de avance: SARS-CoV-2 Mpro apo, NaCl 0.15 M, 310 K",
    )
    story = []

    # Portada: el título va en el lienzo; este spacer deja el bloque marino libre.
    story.append(Spacer(1, 96 * mm))
    story.append(
        P(
            "Enfoque: protocolo de referencia en agua TIP3P con **NaCl 0.15 M** a **310 K**, "
            "modelos de la díada catalítica (H41A, C145A y doble) y despliegue en LARCAD.",
            st["cover_small"],
        )
    )
    story.append(Spacer(1, 6 * mm))
    cover_meta = table(
        ["Campo", "Dato"],
        [
            ["Sistema", "SARS-CoV-2 Mpro apo, cadena A, 306 residuos"],
            ["Condiciones", "AMBER99SB-ILDN / TIP3P / dodecaedro 1.2 nm / 1 bar"],
            ["Sales y T", "NaCl 0.15 M + neutralización · 310 K"],
            ["Producción de referencia", "10 ns (dt = 2 fs, 5 000 000 pasos)"],
            ["Modelos", "silvestre, H41A, C145A, H41A+C145A"],
            ["Clúster", "LARCAD–UNACH · gpu_rtxA5000 (1× RTX A5000) · GROMACS 2026.2"],
            ["Repositorio", "https://github.com/mixcoha/cha-MD-ba"],
            ["Rama / PR", "cursor/6m03-mutantes-h41a-c145a-3810 · PR #14"],
            ["Fecha", "18 de septiembre de 2026"],
            ["Autor", "Edgar Mixcoha"],
        ],
        [4.6 * cm, 12.0 * cm],
        st,
    )
    story.append(cover_meta)
    story.append(PageBreak())

    # --- 1. Proyecto ---
    story.append(P("1. Qué es CHA-MD-BA", st["h1"]))
    story.append(
        P(
            "CHA-MD-BA es un pipeline para preparar, equilibrar, producir y analizar "
            "simulaciones de dinámica molecular de proteínas sobre GROMACS. Hay una "
            "implementación histórica en bash y una versión Python (`cha_md_ba`) con CLI, "
            "API y pruebas. El desarrollo reciente no ha sido añadir módulos genéricos: "
            "ha sido **fijar un benchmark reproducible** sobre un sistema biológico concreto "
            "y correrlo en las mismas condiciones en local y en el clúster LARCAD (UNACH).",
            st["body"],
        )
    )
    story.append(
        bullets(
            [
                "Preparación: PDB → `pdb2gmx` → caja → solvente → iones → minimización.",
                "Equilibración NVT con POSRES decrecientes y NPT a 1 bar.",
                "Producción y análisis (RMSD, RMSF, etc.) en la versión Python.",
                "Artefactos pesados (`.xtc`, `.trr`, `.gro`) en `work/`, fuera de GitHub.",
            ],
            st,
        )
    )

    # --- 2. Por qué este benchmark ---
    story.append(P("2. Por qué este benchmark (6M03)", st["h1"]))
    story.append(
        P(
            "La proteasa principal de SARS-CoV-2 (Mpro, 3CLpro) es la referencia natural "
            "del pipeline: es un blanco terapéutico clásico, tiene estructura apo de alta "
            "calidad (**PDB 6M03**, cadena A, 306 residuos) y una **díada catalítica His41 / Cys145** "
            "cuyo papel se puede contrastar con mutantes a alanina. El benchmark no es un "
            "tutorial vacío: fija condiciones «fisiológicas de laboratorio» — agua explícita, "
            "sal fisiológica y 310 K — para que toda corrida posterior (mutantes, tiempos más "
            "largos, otros campos de fuerza) se compare contra el mismo protocolo.",
            st["body"],
        )
    )
    story.append(
        P(
            "La carpeta del protocolo se llama aún `benchmarks/6M03_nacl_0.5M_310K/` por "
            "historia. **La concentración que usa el código es 0.15 M**, no 0.5 M. El YAML "
            "no se leía solo; la fuente de verdad es `Benchmark6M03Config.ion_concentration = 0.15` "
            "y `gmx genion -neutral -conc 0.15`.",
            st["body"],
        )
    )

    story.append(P("2.1 Pregunta científica que organiza los modelos", st["h2"]))
    story.append(
        P(
            "His41 y Cys145 forman la díada. Truncar cada residuo a alanina (se conservan "
            "N, CA, C, O y CB; `pdb2gmx -ignh` completa hidrógenos) produce tres contrastes "
            "sobre el mismo solvente, la misma caja y el mismo termostato:",
            st["body"],
        )
    )
    story.append(
        table(
            ["Modelo", "Mutación", "run_id", "Para qué"],
            [
                [
                    "Silvestre",
                    "ninguna",
                    "`6M03`",
                    "Referencia: Mpro apo intacta en NaCl 0.15 M / 310 K.",
                ],
                [
                    "1",
                    "H41A",
                    "`6M03_H41A`",
                    "Quita la histidina general de la díada; Cys145 permanece.",
                ],
                [
                    "2",
                    "C145A",
                    "`6M03_C145A`",
                    "Quita el nucleófilo; His41 permanece.",
                ],
                [
                    "3",
                    "H41A + C145A",
                    "`6M03_H41A_C145A`",
                    "Díada completa a alanina; control negativo catalítico.",
                ],
            ],
            [2.4 * cm, 3.0 * cm, 4.2 * cm, 7.8 * cm],
            st,
        )
    )
    story.append(P("Tabla 1. Modelos del benchmark. Numeración PDB, cadena A.", st["caption"]))

    # --- 3. Protocolo ---
    story.append(P("3. Protocolo de referencia", st["h1"]))
    story.append(
        P(
            "Las condiciones están congeladas en `python_version/cha_md_ba/benchmark.py` "
            "(`Benchmark6M03Config`), en `benchmarks/6M03_nacl_0.5M_310K/config.yaml` y en "
            "los `.mdp` de `cluster/larcad/mdp/`. Las pruebas en `test_benchmark_6m03.py` "
            "exigen 310 K, 0.15 M, AMBER99SB-ILDN, TIP3P y `genion -conc 0.15`.",
            st["body"],
        )
    )
    story.append(
        table(
            ["Parámetro", "Valor", "Notas"],
            [
                ["Estructura", "6M03, cadena A, 306 residuos", "Apo; se descartan HOH cristalográficas."],
                ["Campo de fuerzas", "AMBER99SB-ILDN", "`pdb2gmx -ff amber99sb-ildn`"],
                ["Agua", "TIP3P", "Caja dodecaédrica, `-d 1.2` nm al borde."],
                ["Iones", "NaCl 0.15 M + neutro", "`genion -pname NA -nname CL -neutral -conc 0.15`"],
                ["Temperatura", "310 K", "V-rescale; grupos Protein / Non-Protein; τ = 0.1 ps."],
                ["Presión", "1 bar", "NPT de equilibración: C-rescale. Producción: Parrinello–Rahman."],
                ["Paso de tiempo", "2 fs", "`dt = 0.002` ps; LINCS, `constraints = h-bonds`."],
                ["POSRES NVT", "1000 → 800 → 600 → 400 → 200", "kJ mol⁻¹ nm⁻²; 100 ps por constante."],
                ["NPT", "100 ps", "POSRES 200 en el paquete LARCAD; continuación desde NVT."],
                ["Producción", "10 ns", "`nsteps = 5 000 000` en el `.mdp` de referencia."],
                ["HIS", "HIE", "Evita el menú interactivo de protonación en `pdb2gmx`."],
            ],
            [3.4 * cm, 5.4 * cm, 8.6 * cm],
            st,
        )
    )
    story.append(P("Tabla 2. Condiciones del benchmark. Fuente: Benchmark6M03Config y .mdp.", st["caption"]))

    story.append(P("3.1 Tiempos de cada etapa", st["h2"]))
    story.append(
        table(
            ["Etapa", "Integrador", "Pasos", "Tiempo físico", "Salida típica"],
            [
                ["Minimización", "steep", "hasta 50 000", "—", "`minimized.gro`"],
                ["NVT × 5 POSRES", "md / V-rescale", "50 000 c/u", "5 × 100 ps", "`nvt.gro` + `.cpt`"],
                ["NPT", "md / C-rescale", "50 000", "100 ps", "`npt.gro` + `.cpt`"],
                ["Producción", "md / PR", "5 000 000", "10 ns", "`md.xtc` / `md.gro`"],
            ],
            [3.6 * cm, 3.4 * cm, 3.2 * cm, 3.0 * cm, 4.2 * cm],
            st,
        )
    )
    story.append(P("Tabla 3. Duración de referencia. Ampliar producción no cambia el resto del protocolo.", st["caption"]))

    story.append(P("3.2 Composición esperada del sistema", st["h2"]))
    story.append(
        P(
            "Con GROMACS 2023.3, AMBER99SB-ILDN y TIP3P, una corrida de referencia del "
            "silvestre produjo ~21 000 aguas en una caja dodecaédrica de ~700 nm³. Los "
            "conteos de Na⁺/Cl⁻ del YAML histórico (`NA: 215`, `CL: 211`) corresponden a "
            "**0.5 M**; a 0.15 M el número de iones baja y debe leerse de `topol.top` de "
            "esa corrida, no copiarse del YAML.",
            st["body"],
        )
    )
    story.append(
        table(
            ["Cantidad", "Valor de referencia"],
            [
                ["Residuos (cadena A)", "306"],
                ["Átomos de proteína en el PDB (sin H)", "2 367"],
                ["Átomos de proteína con hidrógenos", "4 682"],
                ["Agua TIP3P (SOL)", "~21 061"],
                ["Carga neta tras genion", "0"],
                ["Volumen de caja (dodecaedro)", "~700 nm³"],
            ],
            [8.0 * cm, 8.5 * cm],
            st,
        )
    )
    story.append(P("Tabla 4. Sistema silvestre. Los mutantes a ALA cambian poco el recuento atómico.", st["caption"]))

    # --- 4. Pipeline ---
    story.append(P("4. Pipeline, de PDB a 10 ns", st["h1"]))
    story.append(
        bullets(
            [
                "Descarga de 6M03 desde RCSB (caché local en `work/data` o `cluster/larcad/.cache`).",
                "Limpieza: solo ATOM de proteína; HIS → HIE; se tiran aguas y heteroátomos.",
                "Mutación a ALA si el modelo no es silvestre (`mutate_to_alanine`).",
                "`pdb2gmx` AMBER99SB-ILDN / TIP3P / `-ignh` (y `posre.itp`).",
                "`editconf` dodecaedro centrado, `-d 1.2`.",
                "`solvate` + `grompp` de iones + `genion -neutral -conc 0.15`.",
                "Minimización steepest descent, `emtol = 1000 kJ mol⁻¹ nm⁻¹`.",
                "NVT 310 K, POSRES 1000→200 kJ mol⁻¹ nm⁻² (100 ps cada una).",
                "NPT 310 K / 1 bar / 100 ps (C-rescale + POSRES 200 en LARCAD).",
                "Producción 10 ns, Parrinello–Rahman, sin POSRES.",
            ],
            st,
        )
    )
    story.append(
        P(
            "En local el orquestador es `python scripts/run_benchmark_6m03.py`. En LARCAD "
            "el mismo protocolo está reescrito en `cluster/larcad/run_model.sh` para no "
            "depender de Python en el nodo.",
            st["body"],
        )
    )

    story.append(P("4.1 Cómo lanzarlo en local", st["h2"]))
    story.append(
        Preformatted(
            "python scripts/run_benchmark_6m03.py --resume\n"
            "python scripts/run_benchmark_6m03.py --model 1 --resume    # H41A\n"
            "python scripts/run_benchmark_6m03.py --model 2 --resume    # C145A\n"
            "python scripts/run_benchmark_6m03.py --model 3 --resume    # H41A + C145A\n"
            "python scripts/run_benchmark_6m03.py --model all --resume  # los tres mutantes\n"
            "python scripts/run_benchmark_6m03.py --stages mdps         # solo .mdp a 310 K",
            st["code"],
        )
    )
    story.append(
        P(
            "`--gpu-ids auto` usa la GPU 0 si hay NVIDIA; `--gpu-ids none` fuerza CPU. "
            "La salida vive en `work/<run_id>/` (gitignored).",
            st["body"],
        )
    )

    # --- 5. LARCAD ---
    story.append(P("5. Despliegue en LARCAD", st["h1"]))
    story.append(
        P(
            "El Laboratorio Regional de Cómputo de Alto Desempeño (UNACH) es el destino "
            "de las corridas de 10 ns. Login `ssh -p 212 mixcoha@148.222.27.130`. El "
            "paquete a copiar está en `cluster/larcad/`; se sube con `bash scripts/upload_larcad.sh` "
            "(rsync, no pisa `env.sh` ni `runs/`).",
            st["body"],
        )
    )
    story.append(
        table(
            ["Recurso", "Detalle"],
            [
                ["Cola GPU", "`gpu_rtxA5000` · nodo `gpu1` · 1× RTX A5000"],
                ["Cola CPU", "`larcad` · nodos `n[1-4]`"],
                ["GROMACS GPU", "`gcc-15.1.1` → `mpich-cuda-4.2.3` → `gromacs-mpi-cuda-2026.2`"],
                ["GROMACS CPU", "`gcc-15.1.1` → `mpich-4.2.3` → `gromacs-mpi-2026.2`"],
                ["Binario", "`gmx_mpi`"],
                ["Encolar un modelo", "`bash submit_one.sh 6M03_H41A`"],
                ["Diagnóstico", "`bash diagnose_jobs.sh`"],
            ],
            [4.2 * cm, 12.3 * cm],
            st,
        )
    )
    story.append(P("Tabla 5. Cola y módulos reales de LARCAD (sinfo / module avail).", st["caption"]))

    story.append(
        KeepTogether(
            [
                callout(
                    "Cómo encolar (un modelo; hay una sola GPU)",
                    "En el nodo: `cd ~/cha-md-ba-6m03 && mkdir -p logs && bash submit_one.sh 6M03_H41A`. "
                    "No usar `--export=ALL` ni `sed` sobre `submit_model.slurm`. Los `.out` van al "
                    "directorio de envío (`6M03_H41A_<jobid>.out`). Con una GPU, C145A y el doble "
                    "deben esperar a que H41A corra, o pasarse a la cola CPU `larcad`.",
                    st,
                ),
                Spacer(1, 3 * mm),
            ]
        )
    )

    # --- 6. Estado ---
    story.append(P("6. Estado de las corridas", st["h1"]))
    story.append(
        P(
            "Hasta la fecha de este informe **no hay una producción de 10 ns terminada** "
            "ni en local ni en LARCAD. Lo que sí está cerrado es el protocolo, los PDB "
            "mutados, las pruebas y el paquete de clúster. Lo que sigue es completar las "
            "trayectorias.",
            st["body"],
        )
    )
    story.append(
        table(
            ["Dónde", "Qué se hizo", "Resultado"],
            [
                [
                    "Protocolo / git",
                    "Config 310 K, 0.15 M, mutantes, `.mdp`, tests",
                    "**Hecho.** PR #14 (borrador).",
                ],
                [
                    "PDB mutados",
                    "H41A, C145A, doble (truncamiento a ALA + HIE)",
                    "**Hecho.** En `cluster/larcad/models/`.",
                ],
                [
                    "Local / agente previo",
                    "EM + NVT 100 ps de 6M03, GROMACS 2023.3, 4 CPU, sin GPU",
                    "Parcial (~10.6 ns/día). NPT y 10 ns no se corrieron.",
                ],
                [
                    "LARCAD jobs 1003–1008",
                    "Terna duplicada en `gpu_rtxA5000`",
                    "Muertos. No es la sal.",
                ],
                [
                    "LARCAD jobs 1009–1011",
                    "`sbatch` a mano, `--export=ALL`, `sed` al slurm",
                    "PD → desaparecen en segundos (logs/ + entorno).",
                ],
                [
                    "Producción 10 ns",
                    "Silvestre y tres mutantes",
                    "**Pendiente.**",
                ],
            ],
            [3.3 * cm, 6.4 * cm, 6.8 * cm],
            st,
        )
    )
    story.append(P("Tabla 6. Estado real, sin adelantar trayectorias que no existen.", st["caption"]))

    story.append(P("6.1 Qué no era (el cálculo de la sal)", st["h2"]))
    story.append(
        P(
            "Los jobs no mueren por NaCl 0.15 M. `run_model.sh` llama `genion -conc 0.15`. "
            "Las muertes instantáneas (cola vacía, sin `.err`) coinciden con fallos de "
            "**arranque de SLURM**, no con LINCS a mitad del NPT:",
            st["body"],
        )
    )
    story.append(
        bullets(
            [
                "`#SBATCH --output=logs/%x_%j.out` si `logs/` no existe: Slurm acepta el job y lo mata al abrir el archivo.",
                "`--export=ALL` arrastra MPI/módulos a medias del login al nodo GPU.",
                "Un `sed` sobre `submit_model.slurm` deja el script inconsistente; `env.sh` ya carga gcc + mpich-cuda + GROMACS.",
                "Si el job sí entra a GROMACS: GPU en steep (`-bonded gpu`) y Parrinello–Rahman en el NPT de 100 ps son las causas clásicas de aborto. El paquete actual minimiza en CPU, usa `-nb gpu -pme gpu` solo en dinámica, y equilibra NPT con C-rescale.",
            ],
            st,
        )
    )

    # --- 7. Entregables ---
    story.append(P("7. Entregables de software (lo que ya está en el repo)", st["h1"]))
    story.append(
        table(
            ["Ruta", "Papel en el benchmark"],
            [
                ["`cha_md_ba/benchmark.py`", "Config canónica: T, [NaCl], mutantes, etapas, `--resume`."],
                ["`scripts/run_benchmark_6m03.py`", "Entrada local del benchmark."],
                ["`tests/test_benchmark_6m03.py`", "22 pruebas: 0.15 M, 310 K, HIE, POSRES, modelos."],
                ["`benchmarks/6M03_…310K/`", "README + YAML + `.mdp` de referencia (el nombre aún dice 0.5 M)."],
                ["`cluster/larcad/run_model.sh`", "Mismo pipeline en bash para el nodo."],
                ["`cluster/larcad/mdp/`", "`em`, `nvt`, `nvt_cont`, `npt` (C-rescale), `md` (PR)."],
                ["`models/*/protein.pdb`", "Tres mutantes en `cluster/larcad/models/`, listos para `pdb2gmx`."],
                ["`submit_one.sh`", "Encola sin `--export=ALL`; crea contexto seguro."],
                ["`scripts/upload_larcad.sh`", "rsync al nodo, puerto 212."],
            ],
            [6.6 * cm, 9.9 * cm],
            st,
        )
    )
    story.append(P("Tabla 7. Piezas versionadas. Las trayectorias no se suben a GitHub.", st["caption"]))

    # --- 8. Próximos ---
    story.append(P("8. Próximos pasos", st["h1"]))
    story.append(
        bullets(
            [
                "En LARCAD: un solo job (`6M03_H41A`) hasta ver `pdb2gmx` y NVT en el `.out`. Luego C145A y el doble.",
                "Si `gpu_rtxA5000` sigue saturada, los tres mutantes en paralelo en `larcad` (CPU, GROMACS MPI 2026.2).",
                "No reutilizar `1_preparation` hecha a 0.5 M: borrar `runs/<modelo>` y regenerar iones a 0.15 M.",
                "Cuando existan 10 ns: RMSD/RMSF de proteína, radio de giro y comparación silvestre vs díada mutada.",
                "Opcional: silvestre (`6M03`) en las mismas condiciones, para anclar los tres mutantes.",
                "Ampliar producción (50–100 ns) sin tocar el resto del protocolo, una vez estable el de 10 ns.",
            ],
            st,
        )
    )

    story.append(P("9. Lectura corta", st["h1"]))
    story.append(
        P(
            "El avance de este ciclo **es el benchmark**: un protocolo único (6M03, 0.15 M, "
            "310 K, AMBER99SB-ILDN, TIP3P), cuatro modelos de la díada catalítica, pruebas "
            "automáticas y un paquete SLURM para LARCAD. Las trayectorias de 10 ns son el "
            "siguiente entregable, no el actual. La sal está bien; lo que tumbó los jobs "
            "1003–1011 fue el arranque en el clúster, y eso ya está documentado y parcheado "
            "en la rama del PR #14.",
            st["body"],
        )
    )
    story.append(Spacer(1, 6 * mm))
    story.append(rule())
    story.append(
        P(
            "Documento generado a partir del repositorio `mixcoha/cha-MD-ba`, rama "
            "`cursor/6m03-mutantes-h41a-c145a-3810`. Regenerar: "
            "`python3 scripts/generate_benchmark_progress_pdf.py`.",
            st["caption"],
        )
    )

    def first_page(canvas, doc_):
        cover_page(canvas, doc_)

    def later(canvas, doc_):
        header_footer(canvas, doc_)

    doc.build(story, onFirstPage=first_page, onLaterPages=later)
    return output


if __name__ == "__main__":
    path = build()
    print(f"PDF generado: {path}")
    print(f"Tamaño: {path.stat().st_size} bytes")
