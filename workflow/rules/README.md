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

Las reglas consumidoras combinan el manifest como input `ancient` con su SHA-256
como parámetro. Un mero cambio de timestamp no invalida los análisis. Los
resúmenes y agregados también resuelven directamente el checkpoint: al
reconstruir el grafo se descartan razones de actualización de productores que
ya no deben ejecutarse. Un cambio real del manifest invalida la validación y
los análisis dependientes.

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

`pilot_dada2.smk` consume el QC y los FASTQ recortados, ejecuta inferencia dentro
del estudio y valida las tablas y secuencias de forma independiente. Su
parámetro `config_sha256` invalida DADA2 cuando cambia cualquier byte de la
configuración: coincide con el contrato estricto del registro de hashes de R.
El target `pilot_dada2` termina en `validation.json`, con checksums de salidas y
entradas originales comprobados. El target principal `all` incluye esta etapa.

Las etapas posteriores pueden utilizar `pilot_manifest_path(wildcards)`,
`pilot_manifest_rows(wildcards)`, `pilot_study_id(wildcards)`,
`pilot_run_ids(wildcards)` y `pilot_fastq_paths(wildcards, stage="trimmed")` como
funciones de entrada. Estas funciones consultan el checkpoint, por lo que las
reglas downstream también pueden partir de un repositorio sin datos crudos.

Los tests de integración en `tests/test_workflow.py` comprueban dry-runs en
copias temporales sin manifest, sin FASTQ y con un crudo protegido. Otro test
valida FASTQ sintéticos exclusivamente en un repositorio temporal y comprueba
que modificar sólo el timestamp del manifest no repite trabajo, mientras que
cambiar su contenido sí lo invalida. Ninguna prueba descarga datos ni presenta
los fixtures como evidencia científica.

`pilot_taxonomy.smk` consume DADA2 validado, reconstruye la referencia SILVA
fijada y ejecuta clasificación y validación independiente. El target
`pilot_taxonomy` y `all` requieren `results/taxonomy/pilot/validation.json`.
Sus productores y consumidores resuelven directamente el checkpoint del
manifest. Un SHA-256 de la configuración taxonómica obliga a recalcular
si cambian los parámetros, incluso con timestamp conservado. La configuración
taxonómica separada conserva las entradas del DADA2 original.
