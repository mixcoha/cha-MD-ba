---
cursor:
  subagentId: "bc-a9955558-c07a-51f8-9058-bf66ba110377"
---

# Redock ciego de 2FYT con AutoDock Vina

Sí se reprodujo la pose cristalográfica. Con la misma caja de toda la proteína, exhaustividad **1024** y semilla 42, la pose de mejor afinidad queda a **0,22 Å** del SAH cristalino y puntúa **−9,251 kcal/mol**.

El barrido corto (exhaustividad 8 a 64) no bastó: el mejor RMSD de ese tramo fue 9,30 Å y la mejor pose rank-1 quedó a 22,47 Å. La función de puntuación sí reconoce el sitio (la pose cristalina da −8,834 kcal/mol y una minimización local, −9,048 kcal/mol a 0,13 Å). Lo que faltaba era muestreo. En Vina 1.2.7 un `max_evals` mayor que cero corta la búsqueda; el esfuerzo se sube con la exhaustividad. La caja no se redujo.

## Ligando

| Campo | Valor |
| --- | --- |
| PDB | 2FYT, resolución 2,00 Å |
| Proteína | PRMT3 humana, cadena A, residuos 238–548 |
| Ligando | S-adenosil-L-homocisteína |
| Nombre de residuo | SAH |
| Cadena y número | A 549 |
| Copias | Una sola. La unidad biológica del autor (BIOMOLECULE 1) es monomérica y aplica solo a la cadena A. No hay una segunda copia del ligando. |
| Fórmula | C14 H20 N6 O5 S (26 átomos pesados) |
| Otros HET | Solo 228 aguas (HOH). No hay iones. Se eliminaron todas las aguas. |

## Preparación

Programas: AutoDock Vina **1.2.7** (binario oficial de CCSB-Scripps, función de puntuación `vina`, no smina ni qvina), Meeko 0.8.0, RDKit 2026.09.2, Open Babel 3.1.1, PROPKA 3.5.1, spyrmsd 0.9.0.

**Receptor.** Se conservaron los registros `ATOM` de la cadena A con altloc en blanco o A. De THR A 324, que tiene dos conformeros a ocupación 0,50, se usó el altloc A y se descartó el B. No hay puentes disulfuro ni huecos entre 239 y 548.

PHE A 238 solo tiene en el cristal N, CA, C y O: la cadena lateral no está resuelta. Meeko no puede asignarle la plantilla de fenilalanina. Se omitió ese residuo. Su átomo más cercano a SAH está a 4,82 Å. El resto de la cadena, SER 239–GLN 548, se preparó con `mk_prepare_receptor.py`.

La protonación sigue el pH 7,4 según PROPKA 3.5 sobre el PDB de la cadena A, con las plantillas de Meeko:

- ASP y GLU desprotonados. El pKa más alto de un ASP es 6,98 (ASP 278, a 4,08 Å de SAH); a pH 7,4 la especie mayoritaria sigue desprotonada y se dejó como ASP. El GLU más alto es GLU 346, pKa 5,63.
- LYS y ARG protonados. CYS como tiol neutro (pKa mínimo 8,75; CYS 281, en el sitio, pKa 14,33).
- HIS neutra, no HIP. El pKa más alto es HIS 496, 7,11.
- Tautómeros por geometría de puentes de hidrógeno. HID (protón en ND1): 243, 247, 442 y 496. HIE (protón en NE2): 269, 334, 393, 438, 481 y 523. HIS 247 está a 3,44 Å de SAH: ND1 apunta a carbonilos de cadena principal y NE2 queda a 3,44 Å del O3′ del ligando, así que NE2 se dejó como aceptor.
- Extremo N, SER 239, como NSER (NH3+). Está a 6,9 Å de SAH.
- Extremo C, GLN 548, no trae OXT en el cristal, así que la plantilla de carboxilato terminal no encaja. Se dejó el carbonilo depositado, con el relleno de Meeko. Está a 15,8 Å de SAH.

**Ligando.** Coordenadas y enlaces del SAH cristalográfico (servidor de modelos de RCSB, cadena A, residuo 549). Open Babel `obabel -p 7.4` añadió hidrógenos a pH 7,4 sin mover átomos pesados (desplazamiento máximo 0,000 Å). El resultado es el zwitterión de carga neta 0: amonio y carboxilato, adenina neutra.

```
Nc1ncnc2c1ncn2[C@@H]1O[C@H](CSCC[C@H]([NH3+])C(=O)[O-])[C@@H](O)[C@H]1O
```

`mk_prepare_ligand.py` escribió el PDBQT con cargas de Gasteiger y los hidrógenos no polares fusionados. Quedan 33 átomos (26 pesados y 7 hidrógenos polares) y 11 enlaces rotables (`TORSDOF 11`).

## Caja y parámetros fijos

Una sola caja, la misma en las cuatro corridas. Límites de los 2455 átomos pesados de la proteína (cadena A, altloc en blanco o A, incluidos los cuatro átomos de PHE 238) más 4,0 Å de margen por cada lado.

| | x | y | z |
| --- | ---: | ---: | ---: |
| Mínimo (Å) | −59,325 | −51,943 | −7,770 |
| Máximo (Å) | 5,747 | −3,596 | 51,069 |
| Centro (Å) | −26,789 | −27,7695 | 21,6495 |
| Tamaño (Å) | 73,072 | 56,347 | 66,839 |

Volumen 275 202 Å³. Vina avisó de que supera 27 000 Å³. El SAH cristalográfico queda dentro, a 18,8 Å de la cara más cercana. La caja no se tocó en ninguna corrida. Las cuatro primeras, con 4 CPU, terminaron en menos de dos minutos; por eso el primer barrido no bajó `num_modes`. Después se subió solo la exhaustividad.

El resto quedó fijo: semilla 42, `num_modes` 20, `scoring vina`, `cpu` 4, espaciado 0,375 Å, `energy_range` 3 kcal/mol y `min_rmsd` 1 Å (valores por defecto). Todas las corridas devolvieron 20 modos. La semilla impresa por Vina fue 42.

## Corridas

RMSD de átomos pesados, con simetría y sin superposición, frente al SAH cristalográfico ya protonado. «Mejor pose» es el modo 1 (mejor afinidad). «Mejor RMSD» es el mínimo entre los 20 modos.

| Exhaustividad | Semilla | Afinidad del modo 1 (kcal/mol) | RMSD del modo 1 (Å) | Mejor RMSD (Å) | Modo de ese RMSD | Afinidad de ese modo | Tiempo (s) |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 8 | 42 | −5,383 | 27,986 | 22,217 | 14 | −4,920 | 15 |
| 16 | 42 | −6,052 | 22,469 | 18,745 | 16 | −5,328 | 24 |
| 32 | 42 | −5,935 | 26,542 | 20,337 | 17 | −5,211 | 44 |
| 64 | 42 | −5,935 | 26,542 | 9,295 | 16 | −5,354 | 79 |
| 256 | 42 | −8,409 | 3,768 | 3,768 | 1 | −8,409 | 295 |
| 1024 | 42 | −9,251 | 0,225 | 0,225 | 1 | −9,251 | 1168 |

Centro y tamaño de la caja, en todas las filas: (−26,789, −27,7695, 21,6495) Å y (73,072 × 56,347 × 66,839) Å.

`obrms` de Open Babel, sin la opción `-m`, reproduce esas cifras (22,4685 Å, 26,5417 Å y 9,2953 Å). La opción `-m` superpone las estructuras y no se usó.

## ¿RMSD menor de 2 Å?

Sí, con exhaustividad 1024. El modo 1 queda a **0,22 Å** y es a la vez la pose de mejor afinidad (−9,251 kcal/mol). El modo 3 de esa misma corrida queda a 2,05 Å (−8,382 kcal/mol) y el modo 2 a 3,71 Å (−8,484 kcal/mol).

Exhaustividad 256 ya entra al sitio, pero no cruza el umbral: el modo 1 queda a 3,77 Å y −8,409 kcal/mol. Hasta 64, ninguna pose bajó de 9,30 Å.

El motivo del fallo inicial es el muestreo, no la preparación. Cada corrida independiente de Vina da un número de pasos fijo por el tamaño del ligando (unos 11 rotables); la exhaustividad solo multiplica cuántas de esas corridas se lanzan dentro de la caja. Con 64 no se visita el pozo nativo. Con 1024, sí. La caja siguió siendo la de toda la proteína.

## Ajuste recomendado

Exhaustividad **1024**, semilla **42**, caja de centro (−26,789, −27,7695, 21,6495) Å y tamaño (73,072 × 56,347 × 66,839) Å.

La pose de mejor afinidad reproduce el cristal: **RMSD 0,22 Å** y **afinidad −9,251 kcal/mol**. Tardó 1168 s en 4 CPU.

## RMSD

spyrmsd 0.9.0, `symmrmsd`, solo átomos pesados, sin centrar y sin minimizar. El grafo usa la conectividad, no el orden de enlace, de modo que los dos oxígenos del carboxilato son equivalentes. Comprobaciones: la pose contra sí misma da 0; un desplazamiento de 10 Å da 10 Å (no hay superposición); intercambiar esos dos oxígenos da 0 por simetría y 0,605 Å si se respeta el índice. `obrms` sin `-m` coincide con spyrmsd en las poses citadas.
