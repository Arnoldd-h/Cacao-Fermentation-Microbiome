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

El primer hito se concentra en metadata y no descarga FASTQ de forma masiva. La
vertical slice local contiene únicamente seis runs del piloto.

## Inventario actual

La consulta sistemática documentada recupera 34 BioProjects. Once están
configurados y cribados en detalle; el inventario contiene 1.573 corridas, con
154 incluidas, 101 pendientes y 1.318 excluidas. La comparación reproducible
selecciona `PRJNA492720` como piloto y `PRJNA627078` es el segundo estudio
primario. Otros 13 proyectos permanecen en revisión manual.

El informe y los vacíos de metadata están en
`reports/dataset_inventory.md`.

## Ejecución

Con Python 3.11 o posterior:

```powershell
python scripts/metadata/discover_candidates.py
python scripts/metadata/build_inventory.py
python scripts/metadata/validate_metadata.py
python scripts/metadata/select_pilot_dataset.py
python scripts/metadata/build_pilot_manifest.py
python scripts/metadata/download_pilot_fastq.py
python scripts/metadata/validate_pilot_fastq.py
python -m unittest discover -s tests -v
```

La búsqueda y el inventario consultan únicamente metadata pública de NCBI y
ENA. Los resultados tabulares se escriben de forma atómica; una falla de red o
un cambio de esquema no debe sobrescribir un inventario válido.

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
python/       librería testeable del inventario
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
`CHANGELOG.md`. No existen todavía resultados biológicos ni un meta-análisis.

La selección de licencia y los metadatos completos de citación están pendientes
para no atribuir autores ni términos legales sin confirmación.
