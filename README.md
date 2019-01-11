# cha-MD-ba
Scripts y programas para enviar y analizar simulaciones de dinámica molecular

Antes que nada, abordaremos un poco sobre las Simunlaciones de Dinámica Molecular.

Las simulaciones de dinámica Molecular (SMD) son escenciales para estudiar el comportamiento de sistemas moleculares biológicos en detalle atómico. Obtentemos información no solo estructural sino energética del sistema que estamos estudiando. 

Los pasos necesarios para realizar una SDM son:

  a) Preparación del sistema
    Esto implica tener una configuración inicial, es decir, un conjunto de coordenadas (X,Y,Z) de átomos que comenzaremos a "mover" en condiciones determinadas para obtener la trayectoria
    
  b) Algoritmo de cálculo
    Necesitamos contar con un motor de cálculo que realice la integración de las ecuaciones de movimiento. Existen algunos motores de cálculo libres como: GROMACS o NAMD y otros de pago como AMBER o CHARMM.
    
  c) Computadora.
    Los algoritmos de cálculo se ejecutan en una computadora. El tiempo de ejecución depende de la cantidad de procesadores, la memoria RAM y el sistema operativo de la computadora. Generalmente los sistemas biológicos cuentan con miles de átomos, por eso se necesita una computadora potente para poder integrar numericamente las ecuaciones de movimiento para cada átomo. Entonces, entre más átomos necesitamos más tiempo para resolver las ecuaciones y conocer la trayectoria. Por ello, se utilizan los centros de súpercomputo donde se pueden utilizar, de manera paralela, cientos de procesadores para poder abordar sistemas con muchos átomos.
    
Volviendo a la PREPARACIÓN DEL SISTEMA, existe una serie de consideraciones que debemos tener en cuenta para que un conjunto de coordenadas de átomos pueda tener sentido biológico. De manera general, esas coordenadas se obtienen de manera experimental. Existen técnicas como la difracción de rayos X o la Resonancia Magnética Nuclear que permiten inferir las coordenadas de una proteina o de los ácidos nucleicos y entonces conocer qué forma tienen. Puedes visitar la página del Protein Data Bank, en donde se publican las estructuras resultas de manera experimental, http://www.rcsb.org/pdb/home/home.do
    
Las coordenadas experimentales de un sistema se agrupan en un archivo denominado PDB. Esos archivos cuentan con información, no sólo estructural sino del experimento con el cual se obtuvieron. También tienen información estructural de las moléculas biológicas como: Estructuras secundarias, modificaciones postraduccionales encontradas. ensambles, secuencia, aminoácidos incompletos o regiones faltantes.
    
 El archivo PDB es utilizado como punto de partida, se debe revisar cuidadosamente para estar seguro que no hace falta algun aminoácido o región, que las modificaciones postraduccionales como puentes disulfuro o glicosilaciones son correctas. En el caso de que algo no esté correcto se deben de subsanar los errores o faltantes para tener un sistema completo y con sentido y coherencia biológicos.
    
Existen herramientas computacionales que nos permiten completar átomos faltantes como PDB2PQR. Aunque este algoritmo tiene como finalidad convertir archivos PDB a PQR y obtener los datos de las cargas y radios de los átomos, nos beneficiamos de su algoritmo. Nos permite conocer el estado de protonación de los aminoácidos de una proteína dado que realiza una optimización de la red intramolecular de puentes de hidrógeno. Para ello, puede completar los átomos faltantes de un aminoácido.
    
Si además tenemos un fragmento de varios aminoácidos faltantes, ya no es suficiente con completar las cadenas laterales, se debe realizar un modelado molecular para tener un sistema biológico completo. Para este fin, existen muchos servidores que puedes utilizar como SwissModel o Modeller que te permiten realizar modelados por homología. Generas un archivo PDB, que debe ser evaluado por herramientas que "califican" el modelo.
    
Si ya tenemos el modelo completo, hay información que debemos de conocer de nuestro sistema como si es una proteína soluble, si es necesaria la presencia de un cofactor, si es una proteína transmembranal, la concentración de sales presentes en un entorno biológico, metales necesarios, presencia de azúcares, lípidos, etc. que permitan describir a nuestro sistema lo más similar a un entorno vivo. Así, los resultados que obtengamos serán más confiables y equivalentes a los experimentales. 
    
Lo primero es elegir el campo de fuerzas con el que vamos a describir a nuestro sistema. El campo de fuerzas es un conjunto de valores que definen de qué manera están unidos y cómo interaccionan los átomos entre sí. Contiene los valores de energía de unión (enlaces, ángulos, dihedros, impropios), y de no unión (cargas, radios de van der Waals). Esos valores son usados por el motor de cálculo para resolver las ecuaciones de movimiento. Junto con el campo de fuerzas, se debe elegir el modelo de agua que será utilizado para solvatar, el modelo de agua viene casi implíscito con el campo de fuerza utilizado aunque puede elegirse entre una variedad de modelos de agua que pueden ayudar a explicar procesos específicos como transporte de iones, solvatación, plegamiento de proteínas. Algunas veces es necesario r
    
Como en los experimentos del laboratorio, nuestro soluto (sistema a estudiar) requiere de un tubo de ensayo para ser disuelto y realizar el experimento. En química computacional utilizamos el mismo principio, usamos una caja de solvatación con condiciones periodicas a la frontera. La caja de solvatación puede tener diferentes geometrías y puede ser literalmente una caja cúbica, triclínica, dodecaédrica, prisma hexagonal. Las condiciones periódicas a la frontera implican la continuidad del sistema, las cajas se ensamblan de tal manera que se puebla todo el espacio con cajas donde cada una lleva nuestro soluto disuelto. Las cajas no son aisladas, interaccionan entre sí, por lo que si una molécula (p.ej. agua o iones) salen de una caja y entran exactamente del lado contrario, manteniédose así la misma cantidad de materia. A cada caja se le denomina celda unitaria.
    
Las celdas unitarias corresponderán a qué tipo de sistema estamos estudiando. Por ejemplo una caja cúbica es ideal para moléculas pequeñas. Las dodecaédricas se usan para estudiar proteínas solubles globulares. Sin embargo si nuestro sistema es transmembranal, los primas rectangulares o hexagonales son los más indicados. En el caso de tubos o DNA los prismas hexagonales también son la mejor opción.
    
Así como en el laboratorio, debemos tener en cuenta la temperatura y la presión de nuestros experimentos. Todas las SDM deben tener un factor que describa estas propiedades. Así, existen algoritmos que permiten que la temperatura o la presión puedan ser simuladas. Esos algoritmos permiten modificar o mantener constante la observable, en el caso de la temperatura se llama Termostato y si se trata de la presión se denomina Barostato. En termodinámica, si se realiza un experimento a presión y temperatura constante, se denomina como ensamble canónico. Podemos variar la presión, la temperatura depende de qué experimento queremos realizar. Todos estos 
    
Una vez que tenemos definido el sistema, está completo, solvatado, se ha elegido el campo de fuerzas óptimo, se ha solvatado y se han elegido las condiciones como sales, membrana, contraiones, metales, se hamn elegido las condiciones de simulación, se realiza la SDM.
    
Se obtiene la trayectoria del sistema y luego se analizará evaluando características geométricas como la estructura secundaria de las proteínas, Radio de Giro, Desviación cuadrática média (RMSD). Algunas propiedades como puentes de hidrógeno, análisis de agrupamiento (cluster), Modos normales, Análisis de componentes principales, evaluaciones energéticas, etc. Esos datos nos revelarán el comportamiento de nuestro sistema durante el tiempo de simulación. Nuestro criterio bioquímico nos permitirá darle sentido biológico a esos números permitiéndonos proponer explicaciones del experimento computacional realizado.

En este repositorio encontraremos una serie de scripts que nos permiten automatizar y homogenizar los pasos para realizar las  SDM, para enviar los cálculos al sistema de colas de YOLTLA así como realizar los análisis iniciales de las trayectorias obtenidas.

Síentete libre de aportar, mejorar, interactuar, comentar el software que está aquí.

Edgar Mixcoha
    
    
    
    
    
    
    
   
    
    
