# Resultados

Los outputs actuales se organizan en `qc/`, `tables/`, `dada2/` y `taxonomy/`.
Se conservan tablas pequeñas, diagnósticos y procedencia; FASTQ e intermedios pesados
quedan fuera de Git. Los resultados técnicos del piloto no incluyen todavía
diversidad, inferencia temporal ni meta-análisis.

## Inventario y diseño

`tables/pilot_dataset_selection.tsv` compara los 34 estudios revisados y
conserva PRJNA492720 como piloto. El inventario incluye 4 estudios con 182 runs,
mantiene 7 pendientes y excluye 23. La selección es un resultado metodológico.

`tables/analysis_units.tsv` conserva el mapa run–estudio–lote–tiempo y las
submuestras; `tables/analysis_design.tsv` resume el diseño longitudinal de los
incluidos. Los 12 lotes son las unidades independientes identificadas; los
182 runs no son 182 réplicas biológicas. La procedencia adjunta documenta
configuración, entradas y versión Git.

## QC

`qc/pilot_download_validation.tsv` registra tamaño, MD5, versión de Python y
commit Git para cada FASTQ de la vertical slice. Los archivos de secuenciación
permanecen en `data/raw/` y no se versionan.

`qc/pilot_fastq_validation.tsv` confirma lectura completa del stream gzip y la
estructura FASTQ de cuatro líneas, con conteos de reads, bases y longitudes. No
sustituye los perfiles de calidad de FastQC y MultiQC.

Los HTML y archivos de trabajo de FastQC/MultiQC bajo `qc/pilot/` son
regenerables y permanecen fuera de Git. `qc/pilot/raw_read_quality.tsv` conserva
las métricas crudas esenciales y comprueba que los conteos coincidan con la
validación estructural previa.

`qc/pilot/cutadapt_summary.tsv`, `trimmed_read_quality.tsv` y
`read_quality_comparison.tsv` registran recorte y comparación antes/después.
`primer_detection.tsv` y `primer_detection_trimmed.tsv` conservan la búsqueda
de primers y construcciones configuradas. Los archivos de procedencia enlazan
parámetros y entradas; ninguna tabla demuestra por sí sola ausencia de
contaminación o identidad taxonómica.

## Piloto DADA2

La ejecución integrada sobre seis runs de PRJNA492720 produjo:

| Métrica | Resultado |
|---|---:|
| Pares de entrada | 176.798 |
| Pares después de filtrado | 165.034 |
| Pares ensamblados | 160.649 |
| Pares sin quimeras | 156.561 (88,55 % de entrada) |
| ASVs sin quimeras | 93 |

Las cifras proceden de `dada2/pilot/summary.tsv` y `read_tracking.tsv`. Ninguna
muestra quedó vacía. La corrida del 2026-10-09 UTC (8 de octubre en Colombia)
utilizó `6221445` y pasó validación independiente de 18 entradas y 18 artefactos.
Conteos ASV, secuencias, seguimiento y resumen son idénticos a la primera corrida.
El workflow integrado y su siguiente dry-run terminaron sin trabajo pendiente.

En `dada2/pilot/`, `asv_counts.tsv` contiene abundancias enteras y
`asv_sequences.tsv`/`.fasta` las secuencias. Los IDs están delimitados por estudio
y ejecución; no son claves globales para combinar ASVs de regiones distintas.
`sequence_length_distribution.tsv`, `merge_diagnostics.tsv`, `error_learning.tsv`
y `error_model_F/R.{pdf,svg,png}` permiten revisar ensamblaje, longitudes y
aprendizaje de errores. `warnings.txt` conserva los avisos para su evaluación.

`run_provenance.yaml`, `config_snapshot.yaml`, `input_checksums.tsv`,
`output_checksums.tsv`, `software_versions.tsv` y `session_info.txt` registran
la trazabilidad. `SUCCESS` indica que terminó el ejecutable; `validation.json`
registra la comprobación independiente de sus artefactos. Desde la raíz en WSL:

```bash
python3 scripts/environment/run_in_environment.py python scripts/dada2/validate_outputs.py --run-dir results/dada2/pilot --check-inputs-root .
```

La opción `--check-inputs-root .` comprueba también los hashes de configuración,
código, metadata, QC y FASTQ originales. Una modificación posterior de esas
entradas exige reconstrucción y validación de la etapa afectada.

Los FASTQ filtrados y objetos RDS se guardan en `intermediate/dada2/pilot/`,
fuera de Git. El intento fallido inicial por identificación de pares SRA se
preserva en `dada2/attempts/sra_identifier_failure.yaml`. El método y sus límites
están en [`reports/dada2_methods.md`](../reports/dada2_methods.md).

La retención no demuestra sucesión conservada ni valida parámetros para otros
estudios. Los conjuntos
V4, V3–V4 y full-length se procesarán por separado; el nuevo estudio PacBio
requiere una ruta propia antes de integración taxonómica.

## Taxonomía del piloto

La primera clasificación del 2026-10-09 sobre `d13c9eb`, con checkout limpio al iniciar,
consume las 93 ASVs y 156.561 lecturas originales de seis muestras. Utiliza
SILVA NR99 138.2, R 4.4.3 y DADA2 1.34.0. El validador independiente confirma
21 entradas y 15 artefactos, sin exclusiones ni warnings de clasificación.
Las salidas finales proceden de una repetición sobre `4e7ce24` después de
ajustar la leyenda; las siete tablas científicas son idénticas por SHA-256.
El estado inicial registrado contiene únicamente las salidas previas que
Snakemake retiró antes de regenerarlas; al iniciar no había cambios sin commit
de código o entradas.

En `taxonomy/pilot/`, `taxonomy_unfiltered.tsv` conserva las llamadas sin
máscara; `taxonomy_bootstraps.tsv`, el soporte original; `taxonomy.tsv`, la
asignación primaria con bootstrap 80; y `taxonomy_sensitivity.tsv`, los
umbrales 80/50 del mismo ajuste. `taxonomy_screening.tsv` conserva las marcas
por ASV y sus lecturas originales. `assignment_coverage.tsv` y
`sample_coverage.tsv` describen cobertura global y por muestra, incluyendo
todas las lecturas en sus denominadores.

Con bootstrap 80, 64/93 ASVs tienen etiqueta de género, asociadas a
99.481 lecturas (63,54 %); con 50 son 79/93 y 104.974 lecturas (67,05 %).
Dos ASVs etiquetadas Chloroplast reúnen 32.719 lecturas; tres etiquetadas
Mitochondria reúnen 149. Las cinco conservadas suman 32.868 (20,99 %).
Aunque sus etiquetas Kingdom son Bacteria en SILVA, deben revisarse como
orgánulos antes de diversidad. No se atribuye especie huésped ni contaminación.

Las figuras `assignment_coverage.{pdf,svg,png}` proceden de esas tablas.
`config_snapshot.yaml`, `input_checksums.json`, `provenance.json`,
`software_versions.tsv` y `session_info.txt` documentan parámetros, hashes,
versiones y estado inicial de Git. `SUCCESS` y `validation.json` distinguen
ejecución terminada de validación independiente. El método y las limitaciones
están en [`reports/taxonomy_methods.md`](../reports/taxonomy_methods.md).
