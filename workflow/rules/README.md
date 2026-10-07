# Reglas modulares

`pilot_qc.smk` implementa descarga, validación de integridad, FastQC/MultiQC,
detección y recorte de primers y QC posterior para el piloto. El inventario
permanece en `workflow/Snakefile`.

El checkpoint `build_pilot_manifest` resuelve dinámicamente el estudio y los
runs: no se abre el manifest durante el parseo. Cada FASTQ crudo tiene una
regla productora, salida protegida y validación de bytes/MD5 antes de su
promoción. Los inputs `ancient` de esa regla evitan reemplazar crudos al cambiar
metadata; el reporte agregado comprueba los archivos existentes con el manifest
actual. Las reglas de QC requieren además validación estructural antes de
procesarlos. Cambiar código o parámetros de un paso invalida sus resultados.

La configuración de primers se elige por el BioProject del manifest y verifica
su `study_id`. El método de recorte actual exige `allow_indels=false`,
`discard_untrimmed=false`, `quality_trimming=false` y coincidencia completa del
primer. Otros valores producen un error explícito. Los reportes Cutadapt se
contrastan con los comandos realmente ejecutados, secuencias, parámetros y
rutas; no se aceptan parámetros configurados como prueba de ejecución.

Los resúmenes y detecciones generan archivos `.provenance.json` con fecha UTC,
commit, estado Git modificado, hashes SHA-256 de configuración, código,
entradas y salidas, comando y versiones ejecutadas. Son artefactos pequeños que
acompañan a las tablas TSV; los FASTQ y reportes voluminosos siguen fuera de Git.

Las etapas posteriores pueden utilizar `pilot_manifest_path(wildcards)`,
`pilot_manifest_rows(wildcards)`, `pilot_study_id(wildcards)`,
`pilot_run_ids(wildcards)` y `pilot_fastq_paths(wildcards, stage="trimmed")` como
funciones de entrada. Estas funciones consultan el checkpoint, por lo que las
reglas downstream también pueden partir de un repositorio sin datos crudos.

Los tests de integración en `tests/test_workflow.py` ejecutan únicamente
dry-runs en copias temporales: uno sin manifest y otro sin datos FASTQ. Se
ejecutan dentro del entorno Linux declarado y no descargan datos.
