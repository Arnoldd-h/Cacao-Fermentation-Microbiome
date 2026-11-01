# Auditoría de repetición taxonómica del primer estudio

Fecha local: 2026-10-31. Datos reales: PRJNA492720, 93 muestras y 782 ASVs.

La repetición motivada por corregir el título de cobertura compara una ejecución
sobre `bb2d484` con otra sobre `da43eeb`. Ambas usan SILVA 138.2, la misma
referencia verificada, parámetros, orden de secuencias, semilla R 20260819 y
dos hilos. No se modificó el algoritmo ni se eligió una corrida por sus resultados.

## Observaciones

| Producto | Resultado de la comparación SHA-256 |
|---|---|
| Clasificación primaria bootstrap 80 | Idéntica |
| Clasificación sin máscara | Idéntica |
| Screening primario | Idéntico |
| Versiones de software | Idénticas |
| Soportes bootstrap y sensibilidad | Diferentes |
| Cobertura global y por muestra | Diferente por la sensibilidad |
| Cuatro tablas bacterianas derivadas | Idénticas |
| Siete tablas de diversidad | Idénticas |

La clasificación primaria conserva 475 ASVs y 1.855.013 lecturas asignadas a
género en ambas corridas. A bootstrap 50, el total pasa de 649 ASVs/2.051.014
lecturas a 648 ASVs/2.051.009 lecturas. Es una diferencia **neta** de una ASV
y cinco lecturas; los hashes originales no permiten identificar qué celdas
bootstrap cambiaron. No se afirma que solo cambiara una celda o una secuencia.
La tabla primaria, la separación y la diversidad no cambian en esta comparación.

Los registros de hashes previos y actuales están en
`results/reproducibility/pacheco_montealegre_2020_colombia/taxonomy_caption_comparison.json`.
La última corrida cronológica se conserva con su procedencia y validación;
no se escoge por favorecer una asignación o una interpretación biológica.

## Alcance de la semilla y limitación pendiente

El [código de DADA2 de Bioconductor 3.20](https://github.com/bioc/dada2/blob/RELEASE_3_20/src/taxonomy.cpp)
genera el muestreo bootstrap desde el RNG de R, pero `get_best_genus` inicializa
un generador C++ con `std::random_device` para resolver empates de puntuación.
Ese desempate no queda fijado por `set.seed`. El proyecto también describe
limitaciones de reproducción en [su discusión #1115](https://github.com/benjjneb/dada2/issues/1115).
El mecanismo es compatible con la variación observada; no se identificaron
los empates concretos de esta corrida. Reducir hilos por sí solo no elimina
esa fuente nativa de aleatoriedad.

La validación de integridad, máscaras, fórmulas y hashes de una ejecución
no garantiza que reclasificar produzca idénticos soportes. La coincidencia
primaria de estas dos corridas tampoco garantiza coincidencia en toda repetición.
No se presentan los soportes como plenamente deterministas.

Antes de utilizar asignaciones próximas a umbrales en inferencia o integración,
debe evaluarse una solución reproducible del desempate nativo en un entorno
aislado, registrar cualquier cambio de dependencia/método y repetir sus
validaciones. No se parchea silenciosamente DADA2, no se cambian umbrales
y no se imputan etiquetas para ocultar esta limitación. El hito actual se
limita al procesamiento y descripción del primer estudio, sin inferencia.
