# Diccionario de datos

Los archivos TSV usan UTF-8, una fila de encabezado, campos desconocidos vacíos
y valores decimales con punto. `include` y `analysis_include` aceptan `true`,
`false` o `pending`.

## `studies.tsv`

Contiene una fila por estudio/BioProject candidato. Conserva título y
descripción depositados, referencia publicada cuando fue verificada, diseño,
tecnología, conteos de runs y decisión de screening. `number_samples` cuenta el
subconjunto 16S candidato, no necesariamente todos los BioSamples del proyecto.

## `runs.tsv`

Contiene todos los runs públicos recuperados desde ENA, incluidos WGS, ITS,
controles y runs excluidos. Registra estrategia, selección, layout, instrumento,
conteos, bytes FASTQ estimados, enlaces, checksums y decisión auditable.

`fastq_ftp`, `fastq_md5` y `fastq_bytes` describen los FASTQ distribuidos por
ENA. `submitted_ftp` y `submitted_format` conservan, sin reinterpretación, los
archivos y formatos enviados por los depositantes. Sus nombres pueden aportar
evidencia explícita de marcador, pero no prueban por sí solos ausencia de
procesamiento previo. Un archivo enviado llamado `16S.fastq.gz` no implica dos
lecturas: el layout se conserva desde `library_layout`; en PacBio puede ser
`SINGLE`. La revisión no descarga ninguno de estos archivos.

## `samples.tsv`

Contiene únicamente runs que cumplen el selector 16S de cada candidato incluido
o pendiente. Las columnas mínimas son:

- `study_id`, `sample_id`, `run_accession`, `bioproject`, `biosample`;
- `country`, `location`, `fermentation_hours`,
  `fermentation_duration_hours`, `relative_time`, `fermentation_stage`;
- `season`, `cacao_variety`, `replicate`, `sequencing_platform`, `marker`,
  `region_16s`, `notes`.

Se añaden `fermentation_batch`, `sampling_stratum`, `time_source`,
`analysis_include` y `exclusion_reason` para trazabilidad. Los límites temporales
son: early `<= 0.33`, mid `> 0.33` y `<= 0.66`, late `> 0.66` y `<= 1.0`.

`fermentation_batch` sólo se construye cuando todos los componentes declarados
en `duration_group_fields` están presentes en el identificador publicado. Si
el lote es desconocido queda vacío; `study_id` no sustituye un lote ausente.
Un tiempo absoluto conocido puede conservarse sin duración, tiempo relativo
ni etapa cuando el final de fermentación no se ha verificado. Los textos
anteriores `unknown`/`not reported` se preservan en las notas de configuración,
sin usarlos como valores observados. El año de publicación no se deduce del
año de depósito en SRA/ENA.

## `exclusion_log.tsv`

Registra decisiones a nivel estudio o run con: `entity_type`, `entity_id`,
`study_id`, `bioproject`, `reason`, `metric`, `threshold`, `decision` y `source`.
`pending` representa revisión incompleta, no una exclusión definitiva.

## `discovery_candidates.tsv`

Registra el resultado completo de la búsqueda sistemática configurada en NCBI,
con agregados técnicos derivados de ENA. `preliminary_screen` es sólo triage
automático: `manual_review_*` nunca equivale a inclusión científica. Un proyecto
debe pasar a `config/datasets.yaml` y recibir revisión de paper, marker, tiempo y
diseño antes de incorporarse a `studies.tsv`.

## `results/tables/pilot_dataset_selection.tsv`

Compara todos los estudios configurados y marca un único piloto entre los
estudios primarios elegibles. Registra elegibilidad, ranking, disponibilidad de
publicación y primers, layout, número de corridas candidatas, fermentaciones,
tiempos, etapas y volumen FASTQ. Los criterios y su orden están declarados en
`config/config.yaml`; un estudio `pending` se describe, pero no recibe ranking.

## `pilot_manifest.tsv`

Define la rebanada vertical previa a la descarga: dos muestras de lotes de
fermentación distintos por etapa `early`, `mid` y `late`. Cada fila conserva los
dos enlaces FASTQ paired-end, los MD5 declarados, bytes estimados y la razón de selección.
El archivo es un manifest de entrada; su existencia no implica que los datos ya
hayan sido descargados ni validados por checksum local.

## `results/qc/pilot_download_validation.tsv`

Registra una fila por FASTQ con URL y ruta local, bytes y MD5 esperados y
observados, estado de descarga, fecha UTC, commit Git y versión de Python. Los
estados válidos distinguen descarga nueva, reanudada o reutilizada; todos exigen
coincidencia exacta antes de aceptar el archivo.

## `results/qc/pilot_fastq_validation.tsv`

Registra lectura completa de cada stream gzip y valida la estructura FASTQ. Sus
métricas son número de reads, bases y longitudes mínima/máxima; `gzip_valid` y
`fastq_valid` solo son verdaderos si el archivo completo pudo leerse sin errores
de CRC, truncamiento, encabezado, separador o longitud de calidad.
