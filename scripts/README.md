# Scripts

Los ejecutables actuales cubren entorno, metadata, diseño longitudinal, QC,
DADA2 y taxonomía. Cada script declara sus entradas, salidas y condiciones de fallo en
`--help`. Los comandos siguientes se ejecutan desde la raíz del repositorio en
Linux/WSL.

## Entorno y workflow

`environment/run_in_environment.py` ejecuta un comando en el prefijo existente
de Micromamba para `cacao-microbiome`. Permite `--prefix` o `CACAO_ENV_PREFIX`
cuando la resolución sea ambigua; no instala ni modifica entornos.
`environment/verify_environment.py` comprueba herramientas y paquetes R y
registra sus versiones resueltas.

```bash
python3 scripts/environment/run_in_environment.py python scripts/environment/verify_environment.py
python3 scripts/environment/run_in_environment.py snakemake --snakefile workflow/Snakefile --cores 2 --dry-run
python3 scripts/environment/run_in_environment.py snakemake --snakefile workflow/Snakefile --cores 2
```

El target por defecto alcanza la validación taxonómica del piloto. La ejecución
DADA2 sobre `6221445` conserva su validación; la taxonomía tiene configuración
separada para no modificar sus entradas científicas. La suite aprobó 111 pruebas
Python en Linux y seis comprobaciones taxonómicas en R; el dry-run integrado
queda sin trabajos pendientes.

## Metadata y unidades de análisis

`metadata/discover_candidates.py` consulta NCBI/ENA;
`metadata/build_inventory.py` aplica las decisiones y reglas de
`config/datasets.yaml`; `metadata/validate_metadata.py` comprueba relaciones,
tiempos y elegibilidad. `select_pilot_dataset.py` compara estudios y
`build_pilot_manifest.py` selecciona los seis runs del piloto.

`metadata/build_analysis_units.py` exporta `analysis_units.tsv`,
`analysis_design.tsv` y procedencia. Conserva estudio, lote, tiempo, estrato y
extracciones: los tiempos y las submuestras no se convierten en réplicas
biológicas independientes. Los pendientes quedan fuera del diseño primario.

```bash
python3 scripts/environment/run_in_environment.py python scripts/metadata/discover_candidates.py
python3 scripts/environment/run_in_environment.py python scripts/metadata/build_inventory.py
python3 scripts/environment/run_in_environment.py python scripts/metadata/validate_metadata.py
python3 scripts/environment/run_in_environment.py python scripts/metadata/build_analysis_units.py
```

Estos comandos trabajan con metadata. El inventario actual comprende
34 estudios: 4 incluidos con 182 runs, 7 pendientes y 23 excluidos.

`metadata/download_pilot_fastq.py` es reanudable mediante archivos `.part`, no
sobrescribe un FASTQ final que falle validación y exige coincidencia de bytes y
MD5 antes de promover una descarga a archivo final.

`metadata/validate_pilot_fastq.py` lee íntegramente los streams gzip y valida
estructura, reads, bases y longitudes. El workflow usa además
`qc/download_pilot_read.py` para reconstruir cada FASTQ ausente por separado.

## QC

`qc/summarize_fastqc.py` cruza el reporte tabular de MultiQC con la validación
integral de los FASTQ y falla ante muestras, nombres o conteos inconsistentes.
`qc/detect_primers.py` examina orientaciones y construcciones configuradas;
`qc/trim_primers.py` ejecuta Cutadapt y `qc/summarize_cutadapt.py` verifica sus
parámetros y conteos. `qc/summarize_trimmed_fastqc.py` resume el QC posterior y
lo compara con las entradas crudas. Las reglas están en
`workflow/rules/pilot_qc.smk`.

```bash
python3 scripts/environment/run_in_environment.py snakemake --snakefile workflow/Snakefile --cores 2 pilot_post_trim_qc
```

## DADA2

`dada2/run_pilot.R` y `dada2/helpers.R` procesan un único estudio Illumina
paired-end con parámetros de configuración. Exportan seguimiento de pares,
conteos/secuencias ASV, diagnósticos de errores y longitudes, warnings y
procedencia. `--validate-only` comprueba entradas y parámetros sin ejecutar
filtrado. Los argumentos de ejecución individual están en
[`reports/dada2_methods.md`](../reports/dada2_methods.md).

`dada2/validate_outputs.py` comprueba independientemente checksums, identidad de
muestras, conservación de conteos, secuencias y resúmenes. Con
`--check-inputs-root .` verifica además las entradas originales respecto de la
raíz del repositorio; falla si han cambiado desde la ejecución.

```bash
python3 scripts/environment/run_in_environment.py snakemake --snakefile workflow/Snakefile --cores 2 pilot_dada2
python3 scripts/environment/run_in_environment.py python scripts/dada2/validate_outputs.py --run-dir results/dada2/pilot --check-inputs-root .
python3 scripts/environment/run_in_environment.py python -m unittest discover -s tests -v
python3 scripts/environment/run_in_environment.py Rscript tests/dada2/test_helpers.R
python3 scripts/environment/run_in_environment.py Rscript tests/dada2/test_filter.R
```

El piloto técnico produjo 93 ASVs y conservó 156.561/176.798 pares sin quimeras;
esas cifras no evalúan sucesión. Diversidad, integración y meta-análisis siguen
sin implementar. Los runs PacBio incluidos
en el inventario necesitan un procesamiento específico y no son compatibles
con este ejecutable paired-end Illumina.

## Taxonomía

`taxonomy/download_reference.py` descarga de forma reanudable la referencia
fijada en `config/taxonomy.yaml`; valida bytes, MD5, gzip y FASTA, y guarda SHA-256
en un registro pequeño. La base completa se conserva ignorada por Git.
`taxonomy/run_taxonomy.py` verifica DADA2 y la referencia, ejecuta
`taxonomy/assign_taxonomy.R`, comprueba que no cambiaron las entradas y genera
procedencia y marcador de éxito. `taxonomy/validate_taxonomy.py` compara
independientemente soportes, máscaras, marcas y conteos originales.

```bash
python3 scripts/environment/run_in_environment.py snakemake --snakefile workflow/Snakefile --cores 2 pilot_taxonomy
python3 scripts/environment/run_in_environment.py python scripts/taxonomy/validate_taxonomy.py
python3 scripts/environment/run_in_environment.py Rscript tests/taxonomy/test_helpers.R
```

Las ASVs y etiquetas originales se preservan, con llamada primaria bootstrap 80
y sensibilidad 50. No se asignan especies ni se excluyen ASVs automáticamente.
Método, referencia y limitaciones: [`reports/taxonomy_methods.md`](../reports/taxonomy_methods.md).
