# Resultados

Los outputs regenerables se organizan en `qc/`, `tables/`, `models/` y
`figures/`. `tables/pilot_dataset_selection.tsv` es un resultado metodológico
del inventario, no un resultado biológico. No existen resultados biológicos en
este hito.

`qc/pilot_download_validation.tsv` registra tamaño, MD5, versión de Python y
commit Git para cada FASTQ de la vertical slice. Los archivos de secuenciación
permanecen en `data/raw/` y no se versionan.

`qc/pilot_fastq_validation.tsv` confirma lectura completa del stream gzip y la
estructura FASTQ de cuatro líneas, con conteos de reads, bases y longitudes. No
sustituye los perfiles de calidad que generarán FastQC y MultiQC.

Los HTML y archivos de trabajo de FastQC/MultiQC bajo `qc/pilot/` son
regenerables y permanecen fuera de Git. `qc/pilot/raw_read_quality.tsv` conserva
las métricas crudas esenciales y comprueba que los conteos coincidan con la
validación estructural previa.
