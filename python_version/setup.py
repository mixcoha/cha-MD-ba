from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", "r", encoding="utf-8") as fh:
    requirements = fh.read().splitlines()

setup(
    name="cha-md-ba",
    version="0.1.0",
    author="Edgar Mixcoha",
    author_email="your.email@example.com",
    description="Molecular Dynamics Simulation Pipeline",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/cha-md-ba",
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
    install_requires=requirements,
    extras_require={
        "dev": [
            "pytest>=6.0",
            "black>=21.0",
            "isort>=5.0",
            "mypy>=0.9",
            "flake8>=3.9",
        ],
    },
    entry_points={
        "console_scripts": [
            "cha-md-ba=cha_md_ba.cli:main",
        ],
    },
)

def download_pdb(pdb_code: str, output_path: str) -> bool:
    url = f"https://files.rcsb.org/download/{pdb_code.upper()}.pdb"
    try:
        import requests
        r = requests.get(url)
        if r.status_code == 200:
            with open(output_path, "w") as f:
                f.write(r.text)
            return True
        else:
            print(f"❌ No se pudo descargar {pdb_code} desde el PDB. Código HTTP {r.status_code}")
            return False
    except Exception as e:
        print(f"❌ Error al intentar descargar: {e}")
        return False


def main():
    import sys
    import os
    import shutil

    if len(sys.argv) != 3:
        print("Uso: cha-md-prepare <nombre_simulacion> <archivo.pdb o código PDB>")
        return

    sim_name = sys.argv[1]
    pdb_input = sys.argv[2]
    sim_dir = os.path.join(os.getcwd(), sim_name)
    os.makedirs(sim_dir, exist_ok=True)

    if os.path.isfile(pdb_input):
        pdb_path = os.path.join(sim_dir, os.path.basename(pdb_input))
        shutil.copy(pdb_input, pdb_path)
        print(f"✅ Archivo PDB copiado a {pdb_path}")
    else:
        pdb_path = os.path.join(sim_dir, f"{pdb_input}.pdb")
        if download_pdb(pdb_input, pdb_path):
            print(f"✅ Archivo PDB descargado desde RCSB a {pdb_path}")
        else:
            print("❌ No se pudo obtener el archivo PDB. Abortando.")
            return

    print(f"🏗️ Preparación lista para: {sim_name} con {pdb_path}")