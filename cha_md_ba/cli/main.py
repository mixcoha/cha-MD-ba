"""Módulo CLI para la interfaz de línea de comandos."""

import click
from pathlib import Path

@click.group()
def cli():
    """CHA-MD-BA: Herramienta para análisis de dinámica molecular."""
    pass

@cli.command()
@click.argument('pdb_file', type=click.Path(exists=True))
def prepare(pdb_file):
    """Preparar un sistema para simulación MD."""
    click.echo(f"Preparando sistema desde {pdb_file}")

@cli.command()
@click.argument('input_file', type=click.Path(exists=True))
def simulate(input_file):
    """Ejecutar simulación MD."""
    click.echo(f"Ejecutando simulación con {input_file}")

@cli.command()
@click.argument('trajectory_file', type=click.Path(exists=True))
def analyze(trajectory_file):
    """Analizar trayectoria MD."""
    click.echo(f"Analizando trayectoria {trajectory_file}")

if __name__ == '__main__':
    cli() 