"""Configuración del paquete CHA-MD-BA."""

from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = fh.read().splitlines()

setup(
    name="cha-md-ba",
    version="0.1.0",
    author="Tu Nombre",
    author_email="tu.email@ejemplo.com",
    description="Paquete para análisis de dinámica molecular de proteínas",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/tu-usuario/cha-md-ba",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.20.0",
        "torch>=1.9.0",
        "biopython>=1.79",
        "matplotlib>=3.4.0",
        "seaborn>=0.11.0",
        "pandas>=1.3.0",
        "scipy>=1.7.0",
        "scikit-learn>=0.24.0",
        "tqdm>=4.62.0",
        "fair-esm>=2.0.0",
    ],
    extras_require={
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=4.0",
            "black>=21.0",
            "isort>=5.0",
            "mypy>=0.9",
            "flake8>=3.9",
            "pre-commit>=2.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "cha-md-ba=cha_md_ba.cli.main:cli",
        ],
    },
) 