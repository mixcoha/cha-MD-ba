#!/usr/bin/env python3
"""
Script para generar un PDF con el texto del grant y ejemplos de código de CHA-MD-BA
"""

import os
import re
from reportlab.lib.pagesizes import letter
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Preformatted, PageBreak, Table, TableStyle
from reportlab.lib.colors import Color, HexColor
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY, TA_LEFT
from reportlab.lib import colors
from reportlab.pdfbase.pdfmetrics import registerFontFamily
from reportlab.pdfbase.ttfonts import TTFont

# Crear estilos
styles = getSampleStyleSheet()
styles.add(ParagraphStyle(
    name='CustomTitle',
    fontName='Helvetica-Bold',
    fontSize=20,
    spaceAfter=40,
    spaceBefore=20,
    alignment=TA_CENTER,
    leading=24
))
styles.add(ParagraphStyle(
    name='CustomHeading1',
    fontName='Helvetica-Bold',
    fontSize=16,
    spaceAfter=20,
    spaceBefore=20
))
styles.add(ParagraphStyle(
    name='CustomHeading2',
    fontName='Helvetica-Bold',
    fontSize=14,
    spaceAfter=15,
    spaceBefore=15
))
styles.add(ParagraphStyle(
    name='CustomBodyText',
    fontName='Helvetica',
    fontSize=12,
    spaceAfter=12,
    alignment=TA_JUSTIFY
))
styles.add(ParagraphStyle(
    name='CustomCode',
    fontName='Courier',
    fontSize=10,
    spaceAfter=12,
    backColor=HexColor('#f5f5f5'),
    borderPadding=5,
    borderWidth=1,
    borderColor=HexColor('#e0e0e0')
))
styles.add(ParagraphStyle(
    name='CustomListItem',
    fontName='Helvetica',
    fontSize=12,
    spaceAfter=8,
    leftIndent=20,
    alignment=TA_LEFT
))
styles.add(ParagraphStyle(
    name='CustomNumberedItem',
    fontName='Helvetica',
    fontSize=12,
    spaceAfter=8,
    leftIndent=20,
    alignment=TA_LEFT
))

# Función para convertir texto con formato Markdown a formato ReportLab
def convert_markdown_to_reportlab(text):
    # Convertir negritas (**texto**)
    while '**' in text:
        start = text.find('**')
        end = text.find('**', start + 2)
        if end == -1:
            break
        bold_text = text[start+2:end]
        text = text[:start] + f'<b>{bold_text}</b>' + text[end+2:]
    
    return text

# Leer el contenido del grant
with open('grant_description_technical.md', 'r') as f:
    grant_text = f.read()

# Ejemplos de código
code_examples = {
    'CLI': '''# Ejemplo de uso de la interfaz de línea de comandos
cha-md-ba prepare 1tim.pdb output_dir --forcefield amber99sb-ildn --water-model tip3p
cha-md-ba minimize system.gro system.top min_dir
cha-md-ba nvt system.gro system.top nvt_dir
cha-md-ba npt system.gro system.top npt_dir
cha-md-ba production system.gro system.top prod_dir --num-runs 3
cha-md-ba analyze trajectory.xtc topology.tpr --selection "protein"''',
    
    'Preparación': '''# Preparación de un sistema para simulación MD
from cha_md_ba.prepare import MDSystemPreparator

# Inicializar el preparador
preparator = MDSystemPreparator("1tim.pdb", forcefield="amber99sb-ildn", water_model="tip3p")

# Preparar el sistema
files = preparator.prepare_system(
    output_dir="output_dir",
    box_type="dodecahedron",
    box_size=1.0,
    ions=True,
    minimize=True,
    gpu_ids="0"
)

# Acceder a los archivos generados
print(f"Archivo de topología: {files['topology']}")
print(f"Archivo de estructura: {files['structure']}")''',
    
    'Análisis': '''# Análisis de trayectorias de simulación MD
from cha_md_ba.analysis import MDTrajectoryAnalyzer
import matplotlib.pyplot as plt

# Inicializar el analizador
analyzer = MDTrajectoryAnalyzer("trajectory.xtc", "topology.tpr")

# Calcular RMSD
times, rmsd = analyzer.calculate_rmsd(selection="protein", ref_frame=0)

# Visualizar resultados
analyzer.plot_rmsd("rmsd_plot.png")

# Análisis de enlaces de hidrógeno
hbonds = analyzer.analyze_hbonds("protein", "protein")
print(f"Número promedio de enlaces de hidrógeno: {hbonds['average']}")

# Análisis de clusters
clusters = analyzer.cluster_analysis(selection="protein", cutoff=0.2)
print(f"Número de clusters: {clusters['num_clusters']}")'''
}

# Función para convertir Markdown a elementos de ReportLab
def markdown_to_elements(markdown_text):
    elements = []
    
    # Procesar el texto por líneas
    lines = markdown_text.split('\n')
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        # Título principal
        if line.startswith('# '):
            # Dividir el título en múltiples líneas si es necesario
            title_text = line[2:]
            if len(title_text) > 60:  # Si el título es muy largo
                words = title_text.split()
                lines_title = []
                current_line = ""
                for word in words:
                    if len(current_line + " " + word) <= 60:
                        current_line += " " + word if current_line else word
                    else:
                        lines_title.append(current_line)
                        current_line = word
                if current_line:
                    lines_title.append(current_line)
                
                # Añadir cada línea del título
                for title_line in lines_title:
                    elements.append(Paragraph(title_line, styles['CustomTitle']))
            else:
                elements.append(Paragraph(title_text, styles['CustomTitle']))
            i += 1
        
        # Subtítulos
        elif line.startswith('## '):
            elements.append(Paragraph(line[3:], styles['CustomHeading1']))
            i += 1
        
        # Subtítulos de nivel 2
        elif line.startswith('### '):
            elements.append(Paragraph(line[4:], styles['CustomHeading2']))
            i += 1
        
        # Listas numeradas
        elif re.match(r'^\d+\.\s', line):
            # Convertir el formato Markdown a formato ReportLab
            formatted_line = convert_markdown_to_reportlab(line)
            elements.append(Paragraph(formatted_line, styles['CustomNumberedItem']))
            i += 1
        
        # Listas con viñetas
        elif line.startswith('- '):
            # Convertir el formato Markdown a formato ReportLab
            formatted_line = convert_markdown_to_reportlab(line)
            elements.append(Paragraph(formatted_line, styles['CustomListItem']))
            i += 1
        
        # Texto normal
        elif line:
            # Convertir el formato Markdown a formato ReportLab
            formatted_line = convert_markdown_to_reportlab(line)
            elements.append(Paragraph(formatted_line, styles['CustomBodyText']))
            i += 1
        
        # Línea en blanco
        else:
            elements.append(Spacer(1, 12))
            i += 1
    
    return elements

# Generar el PDF
def generate_pdf(output_filename):
    doc = SimpleDocTemplate(
        output_filename,
        pagesize=letter,
        rightMargin=72,
        leftMargin=72,
        topMargin=72,
        bottomMargin=72
    )
    
    # Contenido del documento
    elements = []
    
    # Convertir el texto del grant a elementos
    elements.extend(markdown_to_elements(grant_text))
    
    # Agregar ejemplos de código
    elements.append(PageBreak())
    elements.append(Paragraph("Ejemplos de Código", styles['CustomHeading1']))
    elements.append(Spacer(1, 12))
    
    for title, code in code_examples.items():
        elements.append(Paragraph(f"Ejemplo: {title}", styles['CustomHeading2']))
        elements.append(Preformatted(code, styles['CustomCode']))
        elements.append(Spacer(1, 20))
    
    # Construir el PDF
    doc.build(elements)
    print(f"PDF generado exitosamente: {output_filename}")

if __name__ == "__main__":
    generate_pdf("CHA-MD-BA_Grant_Description.pdf") 