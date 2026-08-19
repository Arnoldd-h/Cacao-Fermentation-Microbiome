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

## `exclusion_log.tsv`

Registra decisiones a nivel estudio o run con: `entity_type`, `entity_id`,
`study_id`, `bioproject`, `reason`, `metric`, `threshold`, `decision` y `source`.
`pending` representa revisión incompleta, no una exclusión definitiva.
