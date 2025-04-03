Contribuir a CHA-MD-BA
====================

Gracias por su interés en contribuir a CHA-MD-BA. Este documento proporciona pautas e instrucciones para contribuir al proyecto.

Código de Conducta
----------------

Al participar en este proyecto, usted acepta cumplir con nuestro Código de Conducta. Por favor, sea respetuoso y considerado con los demás.

Cómo Contribuir
-------------

1. Haga un fork del repositorio
2. Cree una rama de características
3. Realice sus cambios
4. Envíe un pull request

Configuración del Entorno de Desarrollo
------------------------------------

1. Clone su fork:

.. code-block:: bash

   git clone https://github.com/su_usuario/cha-md-ba.git
   cd cha-md-ba

2. Cree un entorno virtual:

.. code-block:: bash

   python -m venv venv
   source venv/bin/activate  # En Windows: venv\Scripts\activate

3. Instale las dependencias de desarrollo:

.. code-block:: bash

   pip install -e ".[dev]"

Estándares de Código
------------------

- Siga la guía de estilo PEP 8
- Use anotaciones de tipo
- Escriba docstrings para todas las funciones y clases públicas
- Mantenga las funciones enfocadas y pequeñas
- Escriba pruebas para nuevas características

Pruebas
-------

Ejecute las pruebas con:

.. code-block:: bash

   pytest

Documentación
------------

- Actualice la documentación al agregar nuevas características
- Mantenga los docstrings actualizados
- Agregue ejemplos para nuevas características
- Actualice el README si es necesario

Proceso de Pull Request
---------------------

1. Asegúrese de que su código pase todas las pruebas
2. Actualice la documentación
3. Agregue pruebas apropiadas
4. Envíe el pull request con descripción de los cambios

Reporte de Problemas
------------------

Al reportar problemas, por favor incluya:
1. Descripción del problema
2. Pasos para reproducir
3. Comportamiento esperado
4. Comportamiento actual
5. Detalles del entorno

Solicitudes de Características
---------------------------

Para solicitudes de características, por favor:
1. Describa la característica
2. Explique el caso de uso
3. Proporcione ejemplos si es posible

Proceso de Lanzamiento
--------------------

1. Actualice el número de versión
2. Actualice el registro de cambios
3. Cree etiqueta de lanzamiento
4. Construya y suba a PyPI

¿Preguntas?
----------

No dude en abrir un issue o contactar a los mantenedores. 