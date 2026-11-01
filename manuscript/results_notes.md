# Notas de resultados

Actualización: 2026-10-31. Son notas de trabajo basadas en outputs observados;
no constituyen todavía un contraste de hipótesis ni un meta-análisis.

## Inventario y unidad analítica

El inventario incluye cuatro estudios, 182 runs y 12 lotes de fermentación.
Los runs y estratos no equivalen a lotes independientes. Las decisiones y
fuentes están en [estado del proyecto](../PROJECT_STATE.md) y
[unidades analíticas](../results/tables/analysis_units.tsv).

## Piloto técnico de PRJNA492720

Se procesaron seis muestras de cuatro lotes y distintos estratos. De 176.798
pares de entrada se conservaron 156.561 pares sin quimeras (88,55 %) y 93 ASVs.
La separación de cinco ASVs etiquetadas como orgánulos retiró 32.868 lecturas;
la tabla derivada conserva seis muestras, 88 ASVs y 123.693 lecturas.
Fuentes: [DADA2](../results/dada2/pilot/summary.tsv) y
[balance por muestra](../results/filtering/pilot/sample_retention.tsv).

La clasificación primaria utiliza SILVA 138.2 y bootstrap 80. Las tablas por
familia/género conservan todos los rangos ancestrales y las etiquetas vacías.
Con todas las lecturas retenidas de cada muestra como denominador, las cuatro
muestras tomadas entre 36 y 96 horas presentan estas proporciones de Acetobacter:

| Muestra | Horas | Lecturas del género / retenidas | Proporción |
|---|---:|---:|---:|
| BactSaI036z2 | 36 | 35.401 / 36.508 | 96,97 % |
| BactSaII048z1 | 48 | 19.643 / 25.324 | 77,57 % |
| BactSaI072z1 | 72 | 21.446 / 21.839 | 98,20 % |
| BactAnI096z2 | 96 | 16.928 / 16.958 | 99,82 % |

Son proporciones de lecturas 16S, no abundancias celulares absolutas. Las
muestras pertenecen a diferentes combinaciones de lote/estrato y no reconstruyen
una trayectoria general. En las dos muestras de 0 horas, el 94,78 % y 99,39 %
de las lecturas retenidas carece de etiqueta de género. No se imputa su identidad.
Fuentes: [composición](../results/composition/pilot/sample_composition.tsv) y
[cobertura](../results/composition/pilot/sample_label_coverage.tsv).

La riqueza observada va de 2 a 46 ASVs y la diversidad efectiva de Shannon,
de 1,013 a 6,072. La profundidad varía de 9.250 a 36.508 lecturas retenidas.
Son descripciones sin rarefacción ni pruebas entre etapas; no estiman riqueza
total ni demuestran un efecto temporal.
Fuente: [diversidad alfa](../results/diversity/pilot/alpha_diversity.tsv).

## Ampliación del primer estudio

El perfil completo contiene las seis muestras del piloto. Ambos perfiles son
ejecuciones solapadas del mismo estudio: no se suman como 99 muestras ni se
tratan como estudios distintos. Los identificadores ASV pertenecen a cada
ejecución; no se equiparan únicamente por su nombre.

El inventario conserva 94 runs; 93 tienen streams paired-end completos,
correspondientes a cuatro lotes y 47 observaciones lote-tiempo. El run sin pares
completos se excluye de esta ruta por formato, antes del análisis, conservando
su inclusión en el inventario. La descarga de 186 archivos está validada por
bytes, MD5, gzip/CRC y estructura FASTQ.

Los informes crudos contienen 2.728.920 pares, con 4.955–79.920 por muestra.
La calidad por base pasa en 184 archivos y tiene avisos en ambos mates de
SRR7899641. El módulo de adaptadores marca esos dos mates como `fail` y R2
de SRR7899691 como `warn`. Estos indicadores no prueban contaminación biológica
ni justifican por sí solos una exclusión.
Fuente: [registro de ejecución](../reports/full_study_execution_2026-10-31.md).

El recorrido completo conserva 2.540.412 pares filtrados, 2.395.174 ensamblados
y 2.254.530 sin quimeras (82,62 % de la entrada), con 782 ASVs. Las 93 muestras
se conservan. Ambos modelos de error convergen en seis rondas. La validación
DADA2 comprueba 192 entradas y 18 artefactos.
Fuente: [resumen completo DADA2](../results/dada2/pacheco_montealegre_2020_colombia/full/summary.tsv).

La clasificación primaria asigna género a 475/782 ASVs, vinculadas a 1.855.013
lecturas antes de la separación. La sensibilidad a bootstrap 50 asigna género
a 648/782 ASVs y 2.051.009 lecturas en la última corrida; no sustituye la
clasificación primaria. La repetición anterior daba 649 ASVs/2.051.014 lecturas.
La tabla primaria y la diversidad resultan idénticas en ambas corridas, pero
los soportes bootstrap no son plenamente reproducibles con la semilla R.
Véase la [auditoría de repetición](../reports/taxonomy_reproducibility_2026-10-31.md).
Se separan 18 ASVs etiquetadas como orgánulos, con 157.550 lecturas.
La tabla derivada conserva 764 ASVs, 2.096.980 lecturas y las 93 muestras.
Fuentes: [cobertura taxonómica](../results/taxonomy/pacheco_montealegre_2020_colombia/full/assignment_coverage.tsv)
y [balance de separación](../results/filtering/pacheco_montealegre_2020_colombia/full/summary.json).

La riqueza observada va de 1 a 75 ASVs, la diversidad efectiva de Shannon de
1 a 30,363 y la profundidad retenida de 986 a 69.751 lecturas. Las primeras dos
componentes principales con pseudoconteo 1 representan el 12,68 % de la
variación CLR. Estos resúmenes no incluyen pruebas temporales ni corrección
por profundidad o diseño repetido; una separación visual no demuestra un
efecto del tiempo.
Fuentes: [diversidad alfa](../results/diversity/pacheco_montealegre_2020_colombia/full/alpha_diversity.tsv)
y [varianza PCA](../results/diversity/pacheco_montealegre_2020_colombia/full/pca_variance.tsv).

## Preguntas aún sin contrastar

- Repetición de trayectorias bacterianas dentro de los lotes y entre estudios.
- Asociación de taxones con el tiempo y las etapas predefinidas.
- Existencia de un núcleo temporal reproducible.
- Efectos temporales y heterogeneidad entre estudios, sin atribuir efectos
  independientes a factores que el diseño no permita separar.
