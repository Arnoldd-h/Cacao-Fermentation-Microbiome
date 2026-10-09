# Inventario de datasets públicos

Revisión final: 2026-10-07 UTC (6 de octubre en Colombia). Las fechas exactas de
recuperación están en los TSV. El historial inicial del inventario se conserva
en Git; este informe describe el estado actualizado.

## Alcance y método

La búsqueda sistemática NCBI BioProject usa el término configurado:

```text
(cacao[All Fields] OR cocoa[All Fields] OR "Theobroma cacao"[Organism])
AND ferment*[All Fields] AND "bioproject sra"[Filter]
```

NCBI aporta los proyectos y ENA todas sus corridas. La elegibilidad se decide
mediante estrategia, marcador, publicación/protocolo, tiempo, tratamiento y
lote. La revisión de metadata no descarga FASTQ. Los nombres de archivos
`submitted_ftp` se conservan como evidencia depositada para separar 16S/ITS
cuando son explícitos; no sustituyen la comprobación de layout ni prueban que
los reads carezcan de procesamiento previo.

## Estado del cribado

Los 34 BioProjects recuperados están configurados y revisados: 4 incluidos,
7 pendientes y 23 excluidos. La cola de 13 candidatos manuales y los cinco
pendientes anteriores recibieron revisión de fuentes; continúan pendientes
sólo los casos con requisitos concretos aún no resueltos. Las diez exclusiones
automáticas también tienen ahora decisión explícita en configuración.

El inventario tiene 2.357 runs: 182 incluidos, 93 pendientes y 2.082 excluidos.
`samples.tsv` contiene 275 candidatos. El registro de decisiones contiene
2.112 filas: 2.105 exclusiones y 7 pendientes a nivel estudio. Una fila del log
no equivale necesariamente a una muestra.

| BioProject incluido | Subconjunto | Runs | Lotes | Tiempos distintos | FASTQ estimado |
|---|---|---:|---:|---:|---:|
| PRJNA492720, Colombia | 16S V4 Illumina paired-end | 94 | 4 | 13 | 377,07 MiB |
| PRJNA627078, México | 16S V3–V4 Illumina paired-end | 60 | 4 | 5 | 1.320,17 MiB |
| PRJEB40850, Costa Rica 2017 | 16S V4 Illumina paired-end, F1/F2 | 12 | 2 | 6 | 105,65 MiB |
| PRJEB57747, Costa Rica 2019 | 16S full-length PacBio SINGLE, F01/F02 | 16 | 2 | 8 | 158,09 MiB |

Los incluidos suman 2.056.243.128 bytes de FASTQ según ENA y 12 lotes. Estos
conteos no equivalen a 182 réplicas independientes. Los cuatro estudios se
procesarán por separado; PacBio requiere una ruta de procesamiento propia.
Los nuevos conjuntos de Costa Rica tienen años/lotes distintos, aunque
comparten laboratorio y sitio, dependencia que deberá evaluarse al integrar.

La [revisión detallada](inventory_review_2026-10-04.md) documenta los 34
proyectos, fuentes primarias, reglas y limitaciones. Hawaii se excluyó al
verificar que sólo se recuperó el día 1 bacteriano. PRJNA1104253 se excluyó del
primario porque el componente natural es WGS y sus amplicones pertenecen a
ensayos controlados SYNCOM/dropout; no se mezcla con fermentación espontánea
natural. Estas decisiones no cambian umbrales temporales ni se basan en
resultados de abundancia.

## Pendientes con requisitos identificados

| BioProject | Runs candidatos | Evidencia que falta |
|---|---:|---|
| PRJEB53853 | 0 | Separación validada 16S/ITS en bibliotecas mixtas |
| PRJEB82327 | 48 | Protocolo/primers y vínculo inequívoco con publicación |
| PRJEB82871 | 16 | Protocolo/primers y vínculo inequívoco con publicación |
| PRJNA420946 | 14 | Reconciliar 120 h del artículo con controles a 144 h y lotes faltantes |
| PRJNA407677 | 4 | Asignación run-tiempo-lote, protocolo y primers |
| PRJNA935329 | 6 | Protocolo original, diseño, lotes, duración y primers |
| PRJNA842267 | 5 | Métodos completos y unidades independientes del control F0 |

Los pendientes no aportan evidencia al primario. En Camerún se mantienen las
observaciones a 144 h sin inventar duración ni eliminar tiempos para hacerlos
coincidir con el artículo. Los campos desconocidos quedan vacíos, incluido lote
cuando no se puede recuperar. Las descripciones originales y las notas de
fuente permanecen auditables.

## Piloto conservado

El ranking recalculado compara los 34 estudios y conserva PRJNA492720 en primer
lugar. Los criterios y su orden están en `config/config.yaml`; la decisión no
está hardcodeada por accession. Tiene 13 tiempos únicos y cuatro eventos de
fermentación. El resto del ranking incluido es Costa Rica 2017, México y Costa
Rica 2019 según los criterios configurados; el layout de los candidatos se
comprueba a nivel run, porque el resumen del BioProject puede ser mixto.

`metadata/pilot_manifest.tsv` conserva seis runs paired-end: dos por etapa,
de lotes distintos dentro de cada etapa. Sus 12 FASTQ suman 25.154.687 bytes.
La revisión verificó que el manifest calculado permanece idéntico. Los archivos
ya descargados localmente coinciden con bytes/MD5 de ENA y pasaron lectura gzip
y validación FASTQ: 353.596 reads, 78.651.325 bases, longitudes 214–223 nt.
Los FASTQ permanecen fuera de Git. Los resultados QC/DADA2 se describen en sus
informes específicos, sin convertir este piloto técnico en una prueba de
hipótesis biológica.

| Evento colombiano completo | Muestras | Intervalo | Tiempos |
|---|---:|---:|---:|
| Antioquia, periodo I | 24 | 0–132 h | 12 |
| Antioquia, periodo II | 26 | 0–144 h | 13 |
| Santander, periodo I | 20 | 0–108 h | 10 |
| Santander, periodo II | 24 | 0–132 h | 12 |

El tiempo relativo colombiano se deriva por evento usando el máximo observado.
Las etapas se calculan desde configuración. En México las tres extracciones por
caja-tiempo son submuestras; en Colombia los estratos son submuestras del evento.
La unidad independiente para inferencia sigue siendo el lote longitudinal.

## Artefactos y validación

- `config/datasets.yaml`: evidencia, fuentes y selectores de los 34 proyectos.
- `metadata/discovery_candidates.tsv`: 34 proyectos, todos configurados.
- `metadata/studies.tsv`: 34 decisiones y atributos de estudio.
- `metadata/runs.tsv`: 2.357 runs, incluyendo campos enviados y decisión.
- `metadata/samples.tsv`: 275 candidatos con ausencias explícitas.
- `metadata/exclusion_log.tsv`: 2.112 decisiones trazables.
- `metadata/data_dictionary.md`: significado y límites de las columnas.
- `results/tables/pilot_dataset_selection.tsv`: comparación para el piloto.
- `metadata/pilot_manifest.tsv`: seis runs de la validación técnica.

Se completaron descubrimiento, reconstrucción y validación independiente. La
suite inicial ejecutó 85 tests: 81 pasaron y 4 integraciones Snakemake
se omitieron en Windows por requerir el entorno Linux declarado. La validación
integrada posterior en Linux aprobó 89 tests, incluidas esas cuatro integraciones
y las comprobaciones independientes de DADA2. El cierre del workflow del 8 de
octubre aprobó 91 tests, con dos regresiones adicionales del grafo. No se
versionaron cachés, secuencias ni textos completos de artículos. La búsqueda
textual reproducible no garantiza sensibilidad bibliográfica absoluta; los
siguientes cierres dependen de las fuentes primarias específicas documentadas
en la revisión detallada.
