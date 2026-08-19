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
| Configurado y revisado | 6 | Decisión detallada en `config/datasets.yaml` |
| Exclusión Fase I: WGS | 7 | Solo shotgun/WGS en las corridas públicas |
| Exclusión Fase I: otra tecnología/objeto | 3 | RNA-seq, transcriptómica u objeto no pertinente |
| Revisión prioritaria | 5 | Amplicón + señal bacteriana/16S + señal temporal |
| Revisión 16S | 6 | Señal 16S, temporalidad aún no verificada |
| Revisión de marcador | 7 | Amplicón, pero marcador bacteriano no verificado |

La cola prioritaria está formada por `PRJEB53853`, `PRJEB82327`,
`PRJEB82871`, `PRJNA627078` y `PRJNA962540`. Esta clasificación es una ayuda
para el cribado, no un resultado de elegibilidad.

## Estudios configurados

| BioProject | Decisión | Runs públicas | Runs candidatas | Motivo principal |
|---|---:|---:|---:|---|
| PRJNA492720 | incluir | 292 | 94 | Fermentación espontánea longitudinal, 16S V4 paired-end y tiempo cada 12 h |
| PRJNA865318 | pendiente | 63 | 21 | Subconjunto bacteriano disponible; temporalidad y publicación incompletas |
| PRJNA1104253 | pendiente | 473 | 16 | Subconjunto verificable controlado in vitro; no elegible aún para el primario on-farm |
| PRJNA552479 | excluir | 14 | 0 | WGS |
| PRJNA1257864 | excluir | 10 | 0 | WGS |
| PRJNA1264670 | excluir | 9 | 0 | WGS |

El inventario detallado contiene 861 corridas: 94 incluidas, 37 pendientes y
730 excluidas. Las 94 corridas del piloto representan 94 muestras del estudio
colombiano, cuatro eventos de fermentación y aproximadamente 395 MB de FASTQ
comprimido según la metadata pública.

## Piloto seleccionado

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
- La consulta es reproducible, pero no garantiza sensibilidad bibliográfica
  absoluta. La próxima iteración debe complementarla con búsqueda por artículos
  y referencias cruzadas de accessions.

## Artefactos reproducibles

- `metadata/discovery_candidates.tsv`: universo y triaje de 34 BioProjects.
- `metadata/studies.tsv`: decisiones y atributos de los seis estudios revisados.
- `metadata/runs.tsv`: 861 corridas con decisión individual.
- `metadata/samples.tsv`: 131 muestras candidatas armonizadas.
- `metadata/exclusion_log.tsv`: 735 exclusiones trazables.

Todos los artefactos se regeneran con los scripts en `scripts/metadata/` y se
comprueban con el validador independiente y los tests unitarios.
