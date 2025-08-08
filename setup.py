"""Configuración del paquete CHA-MD-BA."""

from setuptools import setup, find_packages

setup(
    name="cha-md-ba",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.20.0",
        "pandas>=1.3.0",
        "matplotlib>=3.4.0",
        "MDAnalysis>=2.0.0",
        "click>=8.0.0",
    ],
    entry_points={
        "console_scripts": [
            "cha-md-ba=cha_md_ba.cli.main:cli",
        ],
    },
    author="CHA-MD-BA Team",
    author_email="team@cha-md-ba.org",
    description="Herramienta para análisis de dinámica molecular",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/cha-md-ba/cha-md-ba",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.8",
) 