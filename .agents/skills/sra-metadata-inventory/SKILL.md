---
name: sra-metadata-inventory
description: Reconstruir, cribar y validar inventarios de BioProjects, SRA y ENA para este meta-análisis de fermentación de cacao. Usar cuando se añadan o revisen accessions, se regeneren metadata/studies.tsv, metadata/runs.tsv, metadata/samples.tsv o metadata/exclusion_log.tsv, se evalúe si un estudio público contiene amplicones bacterianos 16S aptos, o se auditen decisiones de inclusión/exclusión antes de descargar FASTQ.
---

# Inventario de metadata SRA/ENA

Construir un inventario reproducible desde metadata pública sin descargar FASTQ
ni inferir campos ausentes.

## Entradas

- `config/config.yaml`: endpoints y límites temporales.
- `config/datasets.yaml`: BioProjects, evidencia manual verificada y reglas de
  selección por dataset.
- Fuentes oficiales NCBI BioProject y ENA Portal API.

## Flujo

1. Revisar `git status`, `protocol/inclusion_criteria.md` y las decisiones
   existentes antes de cambiar clasificaciones.
2. Verificar cada BioProject en NCBI y obtener todos sus runs desde ENA. No usar
   el título como sustituto de estrategia, marker, layout o muestra.
3. Confirmar publicaciones mediante accession explícito o concordancia
   inequívoca. Guardar DOI/PMID/fuente; dejar vacío lo desconocido.
4. Editar `config/datasets.yaml` para añadir evidencia o reglas; no hardcodear
   excepciones silenciosas en los TSV generados.
5. Ejecutar desde la raíz:

   ```bash
   python scripts/metadata/discover_candidates.py
   python scripts/metadata/build_inventory.py
   python scripts/metadata/validate_metadata.py
   python -m unittest discover -s tests -v
   ```

6. Revisar `git diff`, conteos por decisión y almacenamiento FASTQ estimado.
7. Actualizar protocolo, `PROJECT_STATE.md` y changelog cuando cambie una
   clasificación científica relevante.

## Salidas

- `metadata/studies.tsv`: una fila por estudio candidato.
- `metadata/runs.tsv`: todos los runs recuperados, incluso excluidos.
- `metadata/samples.tsv`: subconjuntos 16S incluidos o pendientes.
- `metadata/exclusion_log.tsv`: decisiones explícitas a nivel estudio/run.
- `metadata/discovery_candidates.tsv`: universo reproducible de búsqueda y triage;
  sus prioridades no son inclusiones automáticas.

## Verificaciones obligatorias

- Accessions únicos y con formato válido.
- Relaciones `sample -> run -> study` sin referencias huérfanas.
- `relative_time` igual a horas/duración y dentro de `[0, 1]`.
- Etapas coherentes con `config/config.yaml`.
- Estudios `include` con al menos tres tiempos, FASTQ público y marcador 16S.
- Runs no 16S, controles y tecnologías fuera de alcance registrados en el log.
- Ninguna credencial, FASTQ o archivo grande agregado a Git.

## Criterios de fallo

- Fallar si NCBI/ENA devuelve error tras los reintentos o un candidato no produce
  un registro BioProject válido.
- Fallar si cambió el esquema esperado de ENA, faltan columnas obligatorias o hay
  duplicados/inconsistencias derivadas.
- Marcar `pending`, no `include`, cuando tiempo, marker, primers, diseño o vínculo
  con la publicación no puedan verificarse.
- Nunca descargar FASTQ como parte de este flujo.
