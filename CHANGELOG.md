# Changelog

Los cambios importantes de este proyecto se documentarán en este archivo.

## [Unreleased]

### Added

- Política persistente de control de versiones y trazabilidad científica.
- Estado inicial del proyecto y registro de decisiones metodológicas.
- Reglas de exclusión para datos ómicos pesados, credenciales y artefactos locales.
- Protocolo científico, criterios de inclusión y registro de decisiones.
- Arquitectura reproducible para metadata NCBI/SRA/ENA con configuración externa.
- Skill local `sra-metadata-inventory` y workflow Snakemake inicial.
- Búsqueda sistemática de BioProjects y triaje conservadora de candidatos.
- Inventarios versionados de estudios, muestras, corridas y exclusiones.
- Validador independiente y 17 tests unitarios.

### Changed

- Nombre y descripción inicial del proyecto actualizados a
  **Cacao Fermentation Microbiome**.
- `PRJNA492720` adoptado como piloto con 94 corridas 16S V4 compatibles con la
  publicación; dos estudios permanecen pendientes y tres WGS se excluyen.
- `PRJNA627078` incorporado como segundo estudio primario con 60 corridas 16S
  V3-V4 longitudinales; tres proyectos prioritarios quedan pendientes y
  `PRJNA962540` se excluye por ausencia de una serie temporal 16S.
