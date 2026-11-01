# Cacao Fermentation Microbiome

Proyecto reproducible para evaluar si existe una sucesión bacteriana conservada
durante la fermentación espontánea del cacao entre estudios, regiones y
condiciones tecnológicas diferentes.

## Pregunta científica

¿Existe una trayectoria de sucesión bacteriana conservada durante la
fermentación espontánea del cacao después de considerar la heterogeneidad entre
estudios?

Las hipótesis iniciales proponen una sucesión parcialmente conservada, un
`core temporal microbiome`, un efecto reproducible del tiempo de fermentación y
heterogeneidad geográfica o tecnológica adicional. Los análisis se diseñarán
para intentar falsar estas hipótesis, no para confirmarlas por construcción.

## Alcance

La primera fase incluye exclusivamente amplicones bacterianos del gen 16S rRNA.
Los estudios se procesarán de forma independiente a nivel ASV y se integrarán
posteriormente mediante taxonomía armonizada y meta-análisis de efectos dentro
de estudio. No se combinarán directamente ASVs generadas con primers o regiones
incompatibles.

## Flujo previsto

```text
BioProject/SRA/ENA metadata
        -> inventario y selección auditable
        -> FASTQ por estudio
        -> FastQC / MultiQC / Cutadapt
        -> DADA2 independiente por estudio
        -> taxonomía y armonización
        -> diversidad y modelos temporales dentro de estudio
        -> meta-análisis de efectos aleatorios
        -> heterogeneidad y sensibilidad
```

El inventario consulta metadata pública. La ruta piloto conserva seis runs;
la ruta ampliada audita los 94 de PRJNA492720 y procesa los 93 con pares completos,
con resultados separados y los mismos parámetros científicos.

## Inventario actual

La consulta sistemática recupera 34 BioProjects, todos configurados y revisados:
4 incluidos, 7 pendientes y 23 excluidos. El inventario contiene 2.357 runs,
con 182 incluidos, 93 pendientes y 2.082 excluidos. Los estudios incluidos son
`PRJNA492720`, `PRJNA627078`, `PRJEB40850` y `PRJEB57747`; representan 12 lotes
de fermentación. Las extracciones y estratos de muestreo se conservan como
submuestras, sin contarlos como lotes independientes.

La selección reproducible mantiene `PRJNA492720` como piloto. Sus seis runs
tienen QC y una ejecución integrada DADA2 con 93 ASVs y 156.561 pares
sin quimeras de 176.798 pares de entrada (88,55 %). Se validaron 18 entradas
y 18 artefactos DADA2. La clasificación con SILVA 138.2 pasó validación
independiente de 21 entradas y 15 artefactos: 64/93 ASVs tienen asignación a
género con bootstrap 80 (63,54 % de las lecturas). Se conservan cinco ASVs
marcadas como orgánulos, con 32.868 lecturas, en las tablas originales.
Su separación auditable produce 88 ASVs y 123.693 lecturas para diversidad
descriptiva, conservando las seis muestras. Se exportan índices alfa,
CLR/Aitchison con sensibilidad, Bray–Curtis y PCA, con validación independiente
de fórmulas y geometría. La suite completa más reciente aprueba 161 pruebas Python.
La inferencia
temporal y la integración entre estudios siguen pendientes.

La ampliación a PRJNA492720 tiene un manifest auditado de 93 runs y cuatro
lotes. El inventario conserva 94: uno no ofrece pares completos. El bloqueo
de ENA quedó resuelto el 2026-10-31: el R2 recuperado coincide en bytes y MD5
con el manifest original. No se excluyó ningún run por el fallo de acceso.
El [registro de descarga](results/planning/pacheco_montealegre_2020_colombia/download_status.json)
documenta la disponibilidad y la validación de los archivos locales.
Se conservan 186/186 FASTQ verificados por bytes, MD5, gzip y estructura.
QC y DADA2 del estudio completo están validados: 93 muestras, 782 ASVs y
2.254.530 pares sin quimeras (82,62 % de la entrada). Taxonomía, separación y
diversidad también están validadas: 764 ASVs y 2.096.980 lecturas retenidas,
con todas las muestras. El piloto solapa seis de esas 93 muestras. El
[registro de ejecución](reports/full_study_execution_2026-10-31.md) conserva
las validaciones y los avisos de calidad sin ocultarlos.

Se exportan [composición por familia/género](results/composition/pacheco_montealegre_2020_colombia/full/sample_composition.tsv)
y [figuras temporales por lote/estrato](results/temporal/pacheco_montealegre_2020_colombia/full/temporal_diversity.png).
Las [notas de resultados](manuscript/results_notes.md) describen el aumento de
la fracción etiquetada como Acetobacter entre el primer y último muestreo de
las ocho series, sin atribuir significación o crecimiento absoluto.

La [auditoría de repetición taxonómica](reports/taxonomy_reproducibility_2026-10-31.md)
documenta una limitación: la semilla R no fija todo el desempate nativo de DADA2.
La primaria y la diversidad se repitieron exactamente, pero cambió la sensibilidad
a 50 en una ASV/cinco lecturas netas. La inferencia sigue pendiente.

`PRJEB57747` aporta lecturas PacBio full-length `SINGLE`; requiere una ruta
propia y no se procesa con el piloto Illumina paired-end actual.

Las fuentes y requisitos de los siete pendientes están en el
[informe del inventario](reports/dataset_inventory.md) y la
[revisión detallada](reports/inventory_review_2026-10-04.md). Los parámetros y
límites del piloto se documentan en [métodos DADA2](reports/dada2_methods.md).

## Ejecución

El stack bioinformático se ejecuta en Ubuntu/WSL2 mediante el entorno
`cacao-microbiome`. Las instrucciones de creación, actualización y verificación
están en [`environment/README.md`](environment/README.md).

Desde la raíz del repositorio en Linux/WSL, con el entorno creado, el wrapper
descubre su prefijo registrado en Micromamba. No requiere activación interactiva:

```bash
python3 scripts/environment/run_in_environment.py python scripts/environment/verify_environment.py
python3 scripts/environment/run_in_environment.py snakemake --snakefile workflow/Snakefile --cores 2 --dry-run
python3 scripts/environment/run_in_environment.py snakemake --snakefile workflow/Snakefile --cores 2
```

Para reconstruir el estudio PRJNA492720 completo, con un presupuesto de 4.500 MB:

```bash
python3 scripts/environment/run_in_environment.py snakemake --snakefile workflow/study.Snakefile --cores 8 --resources mem_mb=4500 --rerun-incomplete
```

El alcance y la exclusión por formato están en [métodos del estudio completo](reports/full_study_methods.md).
Se reutilizan crudos y SILVA verificados; los intermedios y resultados ampliados
usan rutas `full`. Los targets `study_qc` y `study_dada2` permiten ejecutar hasta
esas etapas. Un run sin pares completos permanece incluido en el inventario
y queda excluido explícitamente del manifest paired-end.

El workflow reconstruye inventario, diseño longitudinal, selección y manifest;
descarga los FASTQ ausentes del piloto, ejecuta QC, DADA2, taxonomía,
separación bacteriana y diversidad descriptiva, y valida sus salidas.
La base SILVA pesada se reconstruye desde un depósito con checksum.
Para ejecutar el target DADA2 con sus dependencias, comprobar sus artefactos y
entradas, o ejecutar las pruebas:

```bash
python3 scripts/environment/run_in_environment.py snakemake --snakefile workflow/Snakefile --cores 2 pilot_dada2
python3 scripts/environment/run_in_environment.py python scripts/dada2/validate_outputs.py --run-dir results/dada2/pilot --check-inputs-root .
python3 scripts/environment/run_in_environment.py python -m unittest discover -s tests -v
```

La clasificación del piloto y sus parámetros se describen en
[métodos taxonómicos](reports/taxonomy_methods.md):

```bash
python3 scripts/environment/run_in_environment.py snakemake --snakefile workflow/Snakefile --cores 2 pilot_taxonomy
python3 scripts/environment/run_in_environment.py python scripts/taxonomy/validate_taxonomy.py
```

La diversidad descriptiva y la política de ceros están documentadas en
[métodos de diversidad](reports/diversity_methods.md):

```bash
python3 scripts/environment/run_in_environment.py snakemake --snakefile workflow/Snakefile --cores 2 pilot_diversity
python3 scripts/environment/run_in_environment.py python scripts/filtering/validate_bacterial_table.py
python3 scripts/environment/run_in_environment.py python scripts/diversity/validate_diversity.py
```

Las profundidades retenidas van de 9.250 a 36.508 lecturas y la riqueza observada
de 2 a 46 ASVs. Son descripciones sin rarefacción de seis runs de cuatro lotes;
no reconstruyen sus trayectorias ni demuestran diferencias temporales.

Para generar tablas descriptivas por familia/género después de validar la
separación bacteriana, conservando linajes y metadata temporal:

```bash
python3 scripts/environment/run_in_environment.py python scripts/reporting/summarize_composition.py --analysis full_study
python3 scripts/environment/run_in_environment.py python scripts/reporting/summarize_composition.py --analysis full_study --validate-only
```

`--analysis pilot` aplica la misma política al piloto. Los
[métodos de composición](reports/composition_methods.md) describen el denominador,
las etiquetas desconocidas y la validación; estas tablas no contrastan sucesión.

Para visualizar la diversidad y profundidad frente a las horas originales,
conservando lote y estrato, después de validar la diversidad:

```bash
python3 scripts/environment/run_in_environment.py python scripts/reporting/plot_temporal_diversity.py --analysis full_study
python3 scripts/environment/run_in_environment.py python scripts/reporting/plot_temporal_diversity.py --analysis full_study --validate-only
```

Los [métodos temporales descriptivos](reports/temporal_reporting_methods.md)
explican las líneas y los huecos; el perfil `pilot` muestra puntos aislados.

`--check-inputs-root .` comprueba también los hashes de las entradas originales;
un cambio posterior de esas entradas obliga a reconstruir y validar la etapa
afectada. Los comandos individuales de metadata y QC están descritos en
[`scripts/README.md`](scripts/README.md). La búsqueda y el inventario no
descargan secuencias. Las tablas se escriben de forma atómica y conservan su
mtime cuando el contenido no cambia.

## Reproducibilidad

- Configuración externa en `config/`.
- Inventarios tabulares en `metadata/` con fuentes y decisiones explícitas.
- Workflow principal en Snakemake.
- Tests pequeños sin descargas de FASTQ.
- Datos crudos reconstruibles desde accessions, manifests y checksums.
- Decisiones científicas registradas en `protocol/analysis_decisions.md`.

## Estructura

```text
config/       parámetros y candidatos
metadata/     estudios, muestras, runs y exclusiones
protocol/     protocolo, criterios y decisiones científicas
python/       módulos testeables de metadata, diseño y QC
scripts/      entradas ejecutables del workflow
workflow/     reglas de Snakemake
tests/        pruebas unitarias y de humo
data/         datos raw/interim/processed no versionados cuando son pesados
results/      tablas, modelos y figuras regenerables
manuscript/   notas y texto científico en desarrollo
environment/  dependencias reproducibles
```

## Estado

El estado actual está en `PROJECT_STATE.md`; los cambios importantes en
`CHANGELOG.md`. El 2026-10-31 quedó resuelto el bloqueo de ENA: los 186 FASTQ
del manifest de 93 runs están descargados y validados. Su procesamiento,
diversidad descriptiva, composición y figuras temporales están completos;
las salidas conservan configuración y procedencia. El CLR largo completo es
regenerable y queda fuera de Git, junto a los datos pesados.
La integración entre estudios y evaluación de las
hipótesis siguen pendientes; la descripción de un estudio no demuestra
sucesión conservada entre estudios ni ausencia de contaminación. El límite
del RNG taxonómico se registra antes de avanzar a inferencia.

La selección de licencia y los metadatos completos de citación están pendientes
para no atribuir autores ni términos legales sin confirmación.
