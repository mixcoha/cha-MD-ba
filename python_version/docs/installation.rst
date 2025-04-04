Instalación
===========

Requisitos
----------

CHA-MD-BA requiere Python 3.8 o superior. Se recomienda usar un entorno virtual para la instalación.

Instalación desde PyPI
---------------------

La forma más sencilla de instalar CHA-MD-BA es usando pip:

.. code-block:: bash

   pip install cha-md-ba

Instalación desde el código fuente
---------------------------------

Para instalar desde el código fuente, primero clone el repositorio:

.. code-block:: bash

   git clone https://github.com/mixcoha/cha-MD-ba.git
   cd cha-MD-ba

Luego instale las dependencias y el paquete:

.. code-block:: bash

   pip install -r requirements.txt
   pip install -e .

Verificación de la instalación
-----------------------------

Para verificar que la instalación fue exitosa, puede ejecutar:

.. code-block:: bash

   python -c "import cha_md_ba; print(cha_md_ba.__version__)"

Dependencias
-----------

Las dependencias principales incluyen:

* NumPy >= 1.20.0
* MDAnalysis >= 2.0.0
* Rich >= 10.0.0
* Typer >= 0.4.0
* Pathlib >= 1.0.1

Todas las dependencias se instalarán automáticamente al instalar el paquete. 