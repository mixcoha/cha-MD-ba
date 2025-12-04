#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Script de utilidad para limpiar simulaciones de CompuCell3D.

Este script analiza y gestiona las carpetas de simulaciones de CompuCell3D,
permitiendo identificar y eliminar simulaciones incompletas o inactivas
que ocupan espacio en disco.

Funcionalidades:
- Analiza todas las simulaciones en un directorio base
- Identifica simulaciones incompletas (sin archivos VTK o PNG)
- Detecta simulaciones inactivas (sin modificaciones recientes)
- Calcula el espacio ocupado por cada simulación
- Permite eliminar simulaciones incompletas de forma segura

Uso:
    python limpiar_cc3d.py [--dir DIRECTORIO] [--horas HORAS]

Argumentos:
    --dir DIRECTORIO    Directorio base de simulaciones (default: CC3DWorkspace)
    --horas HORAS       Horas de inactividad para considerar una simulación como inactiva (default: 12)
"""

import os
import time
from pathlib import Path
from datetime import datetime
import shutil
import argparse
from typing import List, Dict, Any
import sys

def analizar_simulaciones(directorio_base: str, horas_inactividad: int = 12) -> List[Dict[str, Any]]:
    """
    Analiza todas las simulaciones en el directorio base.
    
    Args:
        directorio_base (str): Ruta al directorio que contiene las simulaciones
        horas_inactividad (int): Horas sin modificaciones para considerar inactiva
        
    Returns:
        list: Lista de diccionarios con información de cada simulación
    """
    ahora = time.time()
    resultados = []
    
    try:
        directorio = Path(directorio_base)
        if not directorio.exists():
            print(f"❌ Error: El directorio {directorio_base} no existe")
            sys.exit(1)
            
        for carpeta in directorio.iterdir():
            if not carpeta.is_dir():
                continue
                
            latticedata = carpeta / "latticedata"
            confield_dirs = [d for d in carpeta.iterdir() if d.is_dir() and d.name.endswith("_COnField")]
            
            vtk_archivos = list(latticedata.glob("*.vtk")) if latticedata.exists() else []
            png_archivos = [png for d in confield_dirs for png in d.glob("*.png")]
            
            ultima_mod = max((f.stat().st_mtime for f in carpeta.rglob("*") if f.is_file()), default=0)
            tiempo_inactivo_horas = (ahora - ultima_mod) / 3600
            tamaño_total = sum(f.stat().st_size for f in carpeta.rglob("*") if f.is_file()) / (1024**3)
            
            if vtk_archivos or png_archivos:
                estado = "✅ Finalizada" if tiempo_inactivo_horas < horas_inactividad else "✅ (Inactiva)"
            else:
                estado = "⚠️ Incompleta" if tiempo_inactivo_horas > horas_inactividad else "⏳ En proceso"
                
            resultados.append({
                "path": carpeta,
                "Simulación": carpeta.name,
                "Estado": estado,
                "VTK archivos": len(vtk_archivos),
                "PNG archivos": len(png_archivos),
                "Última modificación": datetime.fromtimestamp(ultima_mod).strftime("%Y-%m-%d %H:%M:%S"),
                "Inactiva (hrs)": round(tiempo_inactivo_horas, 1),
                "Tamaño (GB)": round(tamaño_total, 2)
            })
            
    except Exception as e:
        print(f"❌ Error al analizar simulaciones: {e}")
        sys.exit(1)
        
    return resultados

def borrar_simulaciones_incompletas(resultados: List[Dict[str, Any]]) -> None:
    """
    Elimina las simulaciones marcadas como incompletas.
    
    Args:
        resultados (list): Lista de resultados del análisis de simulaciones
    """
    incompletas = [r for r in resultados if r["Estado"] == "⚠️ Incompleta"]
    
    if not incompletas:
        print("No se encontraron simulaciones incompletas.")
        return
        
    print("\nSimulaciones incompletas encontradas:")
    for r in incompletas:
        print(f" - {r['Simulación']} ({r['Tamaño (GB)']} GB, inactiva {r['Inactiva (hrs)']}h)")
        
    confirmar = input("\n¿Deseas borrar estas carpetas? (sí/no): ").strip().lower()
    if confirmar != "sí":
        print("Operación cancelada.")
        return
        
    for r in incompletas:
        try:
            shutil.rmtree(r["path"])
            print(f"✔️ Borrado: {r['Simulación']}")
        except Exception as e:
            print(f"❌ Error al borrar {r['Simulación']}: {e}")

def main():
    """Función principal del script."""
    parser = argparse.ArgumentParser(description="Limpieza de simulaciones CompuCell3D")
    parser.add_argument("--dir", default="/Users/mixcoha/CC3DWorkspace/",
                      help="Directorio base de simulaciones")
    parser.add_argument("--horas", type=int, default=12,
                      help="Horas de inactividad para considerar una simulación como inactiva")
    
    args = parser.parse_args()
    
    print(f"Analizando simulaciones en: {args.dir}")
    resultados = analizar_simulaciones(args.dir, args.horas)
    
    print("\nResumen de simulaciones:")
    for r in resultados:
        print(f"{r['Estado']}: {r['Simulación']} ({r['Tamaño (GB)']} GB)")
    
    borrar_simulaciones_incompletas(resultados)

if __name__ == "__main__":
    main() 