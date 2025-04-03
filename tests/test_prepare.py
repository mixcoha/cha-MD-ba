import os
import tempfile
import requests
from unittest.mock import patch, Mock
from cha_md_ba.prepare import download_pdb


@patch("requests.get")
def test_download_pdb_success(mock_get):
    mock_response = Mock()
    mock_response.status_code = 200
    mock_response.text = "HEADER TEST PDB"
    mock_get.return_value = mock_response

    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = os.path.join(tmpdir, "test.pdb")
        success = download_pdb("1tim", output_path)

        assert success
        assert os.path.exists(output_path)
        with open(output_path) as f:
            assert f.read() == "HEADER TEST PDB"

@patch("requests.get")
def test_download_pdb_failure(mock_get):
    mock_response = Mock()
    mock_response.status_code = 404
    mock_get.return_value = mock_response

    with tempfile.TemporaryDirectory() as tmpdir:
        output_path = os.path.join(tmpdir, "test.pdb")
        success = download_pdb("invalidcode", output_path)

        assert not success
        assert not os.path.exists(output_path)

# Crea la carpeta scripts/ si no existe
os.makedirs("scripts", exist_ok=True)

# Dentro de scripts/, crea un archivo llamado run_preparation.py
with open("scripts/run_preparation.py", "w") as f:
    f.write("""import os
from cha_md_ba.prepare import download_pdb, MDSystemPreparator

# Paso 1: Descargar el archivo PDB (si aún no lo tienes)
pdb_id = "1tim"
pdb_path = f"data/{{pdb_id}}.pdb"

# Asegúrate de que el directorio 'data/' existe
os.makedirs("data", exist_ok=True)

# Descargar desde RCSB
if not os.path.exists(pdb_path):
    success = download_pdb(pdb_id, pdb_path)
    if not success:
        raise RuntimeError(f"No se pudo descargar {{pdb_id}} desde RCSB.")
else:
    print(f"Archivo {{pdb_path}} ya existe.")

# Paso 2: Preparar el sistema con GROMACS
preparator = MDSystemPreparator(
    pdb_path=pdb_path,
    forcefield="amber99sb-ildn",
    water_model="tip3p"
)

# Lanzar preparación
result = preparator.prepare_system(output_dir="simulations")

# Verificar resultados
print("\\n🎉 Preparación completa. Archivos generados:")
for key, path in result.items():
    print(f"{{key}}: {{path}}")
# gmx se asume como 'gmx_mpi' por defecto en la clase
""")
