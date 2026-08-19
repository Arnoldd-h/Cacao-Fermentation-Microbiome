# Inventario inicial de datasets públicos

Fecha de recuperación: 2026-08-19 (UTC)

## Alcance y método

Este hito inventaría metadata pública; no descarga archivos FASTQ. La búsqueda
sistemática se ejecutó contra NCBI BioProject con el término documentado:

```text
(cacao[All Fields] OR cocoa[All Fields] OR "Theobroma cacao"[Organism])
AND ferment*[All Fields] AND "bioproject sra"[Filter]
```

NCBI aporta el universo de BioProjects y sus descripciones; ENA aporta la tabla
de corridas, tecnología, layout, tamaño y enlaces de descarga. La triaje
automática es deliberadamente conservadora: nunca convierte una coincidencia en
inclusión. Los proyectos prioritarios deben pasar revisión de publicación,
marcador, diseño temporal y compatibilidad con el protocolo.

## Resultado de la búsqueda sistemática

Se recuperaron 34 BioProjects:

| Clasificación preliminar | Proyectos | Interpretación |
|---|---:|---|
| Configurado y revisado | 11 | Decisión detallada en `config/datasets.yaml` |
| Exclusión Fase I: WGS | 7 | Solo shotgun/WGS en las corridas públicas |
| Exclusión Fase I: otra tecnología/objeto | 3 | RNA-seq, transcriptómica u objeto no pertinente |
| Revisión 16S | 6 | Señal 16S, temporalidad aún no verificada |
| Revisión de marcador | 7 | Amplicón, pero marcador bacteriano no verificado |

Los cinco candidatos que inicialmente formaban la cola prioritaria ya fueron
configurados. La cola sistemática restante contiene 13 proyectos: seis con
señal 16S y siete que requieren comprobar el marcador. Esta clasificación es
una ayuda para el cribado, no un resultado de elegibilidad.

## Estudios configurados

| BioProject | Decisión | Runs públicas | Runs candidatas | Motivo principal |
|---|---:|---:|---:|---|
| PRJNA492720 | incluir | 292 | 94 | Fermentación espontánea longitudinal, 16S V4 paired-end y tiempo cada 12 h |
| PRJNA865318 | pendiente | 63 | 21 | Subconjunto bacteriano disponible; temporalidad y publicación incompletas |
| PRJNA1104253 | pendiente | 473 | 16 | Subconjunto verificable controlado in vitro; no elegible aún para el primario on-farm |
| PRJNA552479 | excluir | 14 | 0 | WGS |
| PRJNA1257864 | excluir | 10 | 0 | WGS |
| PRJNA1264670 | excluir | 9 | 0 | WGS |
| PRJNA627078 | incluir | 132 | 60 | Serie espontánea V3-V4; se excluyen 26S y controles de cosecha limpia |
| PRJEB53853 | pendiente | 60 | 0 | Falta mapear inequívocamente marcador y tratamiento por corrida |
| PRJEB82327 | pendiente | 369 | 48 | Controles espontáneos full-length 16S identificados; paper/primers pendientes |
| PRJEB82871 | pendiente | 140 | 16 | Controles espontáneos full-length 16S identificados; paper/primers pendientes |
| PRJNA962540 | excluir | 11 | 0 | Amplicones sin tiempo y solo dos tiempos WGS |

El inventario detallado contiene 1.573 corridas: 154 incluidas, 101 pendientes y
1.318 excluidas. Las 94 corridas del piloto representan 94 muestras del estudio
colombiano, cuatro eventos de fermentación y aproximadamente 395 MB de FASTQ
comprimido según la metadata pública.

## Selección del piloto

La tabla `results/tables/pilot_dataset_selection.tsv` compara los 11 estudios
mediante criterios declarados en `config/config.yaml`. Los dos estudios primarios
son elegibles, pero `PRJNA492720` ocupa el primer lugar por su serie de 13 tiempos
y menor volumen estimado (377,07 MiB), frente a cinco tiempos y 1.320,17 MiB en
`PRJNA627078`. La decisión no depende de un accession codificado en el script.

`metadata/pilot_manifest.tsv` prepara seis runs paired-end de `PRJNA492720`: dos
por etapa temporal y de fermentaciones distintas dentro de cada etapa. Sus 12
FASTQ suman 25.154.687 bytes estimados. No se han descargado; los MD5 del manifest
son los declarados por ENA y todavía deben verificarse sobre archivos locales.

### Diseño temporal del piloto

`PRJNA492720` es la vertical slice inicial. El subconjunto publicado y
compatible comprende Antioquia y Santander:

| Evento | Muestras | Intervalo | Puntos temporales |
|---|---:|---:|---:|
| Antioquia, periodo I | 24 | 0–132 h | 12 |
| Antioquia, periodo II | 26 | 0–144 h | 13 |
| Santander, periodo I | 20 | 0–108 h | 10 |
| Santander, periodo II | 24 | 0–132 h | 12 |

El tiempo relativo se deriva por evento usando la duración máxima observada y
se valida en el intervalo cerrado `[0, 1]`. Las etapas `early`, `mid` y `late`
se calculan desde configuración; no están codificadas manualmente.

## Vacíos y riesgos de metadata

- Los primers, región 16S, tipo de fermentación y DOI no son uniformemente
  recuperables desde SRA/ENA y requieren verificación en publicaciones.
- Algunos BioProjects mezclan 16S, ITS, WGS, fermentaciones on-farm y ensayos
  inoculados; el filtro debe operar a nivel de corrida y muestra.
- Una búsqueda textual puede producir falsos positivos, por ejemplo estudios
  de levaduras, intestino de ratón o matrices distintas al cacao fermentado.
- `PRJNA492720` contiene más corridas que el subconjunto 16S descrito en el
  artículo; el inventario conserva las exclusiones para evitar selección
  silenciosa.
- En `PRJNA627078`, las tres extracciones por caja y tiempo no deben tratarse
  como réplicas biológicas independientes; la caja es la unidad longitudinal.
- La consulta es reproducible, pero no garantiza sensibilidad bibliográfica
  absoluta. La próxima iteración debe complementarla con búsqueda por artículos
  y referencias cruzadas de accessions.

## Artefactos reproducibles

- `metadata/discovery_candidates.tsv`: universo y triaje de 34 BioProjects.
- `metadata/studies.tsv`: decisiones y atributos de los 11 estudios revisados.
- `metadata/runs.tsv`: 1.573 corridas con decisión individual.
- `metadata/samples.tsv`: 255 muestras candidatas armonizadas.
- `metadata/exclusion_log.tsv`: 1.327 exclusiones trazables.
- `results/tables/pilot_dataset_selection.tsv`: comparación y ranking del piloto.
- `metadata/pilot_manifest.tsv`: seis muestras preparadas para la vertical slice.

Todos los artefactos se regeneran con los scripts en `scripts/metadata/` y se
comprueban con el validador independiente y los tests unitarios.
