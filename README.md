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

El inventario consulta metadata pública. La descarga y el procesamiento local
se limitan a seis runs del piloto; no se descargan todos los proyectos.

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
de fórmulas y geometría. La suite aprueba 125 pruebas Python. La inferencia
temporal y la integración entre estudios siguen pendientes.

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
`CHANGELOG.md`. La integración entre estudios y evaluación de las
hipótesis siguen pendientes; el resultado técnico del piloto no demuestra
sucesión conservada ni ausencia de contaminación.

La selección de licencia y los metadatos completos de citación están pendientes
para no atribuir autores ni términos legales sin confirmación.
