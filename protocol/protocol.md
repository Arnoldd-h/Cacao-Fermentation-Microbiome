# Protocolo científico inicial

## Pregunta y objetivo

Evaluar si existe una sucesión bacteriana conservada durante la fermentación
espontánea del cacao entre estudios, regiones geográficas y condiciones de
fermentación diferentes.

Los objetivos secundarios son identificar taxones asociados a etapas temporales,
definir un `core temporal microbiome` reproducible y cuantificar la variación
atribuible a tiempo, estudio, geografía, temporada, genotipo, sistema de
fermentación y tecnología.

## Hipótesis a falsar

1. Existe una trayectoria de sucesión microbiana parcialmente conservada.
2. Existe un conjunto recurrente de taxones que forma un core temporal.
3. El tiempo explica una fracción reproducible de la estructura comunitaria aun
   considerando estudio y región.
4. La geografía y la tecnología generan heterogeneidad adicional.

Las hipótesis no se modificarán para acomodarlas a los resultados.

## Alcance de la primera fase

- Amplicones bacterianos 16S rRNA con reads crudos públicos.
- Fermentación de granos o masa pulpa-grano de `Theobroma cacao`.
- Diseño longitudinal con tiempo absoluto, relativo o una secuencia temporal
  recuperable.
- Procesamiento independiente por estudio; integración cross-study en género o
  familia según resolución defendible.

ITS y shotgun quedan fuera del análisis primario hasta estabilizar el pipeline
16S. Los datasets mixtos pueden aportar únicamente su subconjunto 16S si la
separación es inequívoca.

## Inventario y selección

Las semillas iniciales son PRJNA492720, PRJNA865318, PRJNA1104253,
PRJNA552479, PRJNA1257864 y PRJNA1264670. No constituyen inclusiones
automáticas.

Para cada proyecto se consultarán NCBI BioProject y ENA, se identificarán runs,
BioSamples, estrategia, selección de librería, plataforma, layout, tamaño
estimado y disponibilidad de FASTQ. Las publicaciones se vincularán solamente
con evidencia explícita del accession o concordancia verificable de diseño y
datos.

El inventario se reconstruye sin descargar FASTQ. Toda decisión se conserva en
`metadata/exclusion_log.tsv` y puede ser `include`, `exclude` o `pending`.

## Unidad analítica y tiempo

`study_id` identifica la publicación o estudio. Cuando un BioProject contiene
múltiples fermentaciones, cada combinación de ubicación/temporada/lote se
conserva como `fermentation_batch`.

La unidad independiente es el lote de fermentación; el tiempo es una medida
repetida dentro del lote. Extracciones y estratos de una misma combinación
estudio/lote/hora se conservan como submuestras, sin contarlas como fermentaciones
independientes. `results/tables/analysis_units.tsv` registra esta jerarquía y
`analysis_design.tsv` resume sus tamaños antes de cualquier modelo. La futura
agregación de abundancias dentro de cada observación se implementará y validará
por estudio; estas tablas no suman ni modifican conteos ASV.

La rebanada técnica de seis runs no se usará para contrastar sucesión, estimar
efectos temporales ni construir un core. En los dos estudios iniciales país,
estudio y región 16S están asociados; no se atribuirá un efecto independiente a
geografía o tecnología sin un diseño que lo identifique. La ampliación del
inventario obliga a reevaluar esa identificabilidad, no la garantiza por sí sola.

Se mantienen `fermentation_hours` y `fermentation_duration_hours`. El tiempo
relativo se calcula únicamente cuando ambos son conocidos y la duración es
positiva:

```text
relative_time = fermentation_hours / fermentation_duration_hours
```

La clasificación inicial usa límites configurables: `early` hasta 0.33
inclusive, `mid` mayor que 0.33 y hasta 0.66 inclusive, y `late` mayor que 0.66.
El análisis conservará además `relative_time` continuo y evaluará sensibilidad a
los límites.

Cuando la duración procede del último tiempo observado, se registra esa fuente
y se contrasta el análisis de tiempo relativo con horas absolutas. No se supone
que el último muestreo representa el mismo estado biológico final en todos los
lotes.

## Procesamiento previsto

1. Inspección de metadata y manifest por estudio.
2. Estimación de almacenamiento antes de descargar.
3. FastQC y MultiQC sobre reads crudos.
4. Inspección y remoción documentada de primers con Cutadapt cuando corresponda.
5. DADA2 independiente con parámetros derivados de la calidad de cada estudio.
6. Asignación taxonómica con una versión registrada de SILVA.
7. Separación de cloroplasto, mitocondria, no bacterias y contaminantes bajo
   reglas explícitas.
8. Exportación por estudio de ASVs, secuencias, taxonomía, QC y metadata.

## Análisis estadístico previsto

- Alfa diversidad dentro de estudio: riqueza observada, Shannon y Simpson.
- Beta diversidad primaria CLR/Aitchison y secundaria Bray-Curtis.
- PCA/PCoA, PERMANOVA y análisis de dispersión antes de interpretar efectos.
- Modelos temporales por estudio con horas, tiempo relativo y etapa.
- Abundancia diferencial composicional, inicialmente ANCOM-BC2 si los datos y
  el diseño lo permiten, con FDR.
- Efectos por taxón dentro de estudio combinados mediante efectos aleatorios
  REML en `metafor`, reportando IC 95 %, Q, tau², I² y número de estudios.
- Sensibilidad leave-one-study-out, límites temporales, core, resolución
  género/familia y tratamiento de taxones raros.

No se afirmará que existe un meta-análisis hasta disponer de varios estudios
independientes apropiados.

## Exclusiones y control de calidad

Ninguna muestra se excluirá silenciosamente. Cada decisión deberá registrar
entidad, razón, métrica, umbral y estado. El QC por muestra conservará lecturas
de entrada, filtradas, denoised, merged, non-chimeric y porcentaje de retención.

## Trazabilidad

Cada resultado científico importante deberá registrar fuentes de datos,
configuración, versiones de software, fecha y commit Git. Resultados observados
e interpretación se mantendrán separados.
