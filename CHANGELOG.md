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
- Validador independiente y 30 tests unitarios.
- Comparación reproducible del piloto y manifest equilibrado de seis runs.
- Descarga reanudable con comprobación de espacio, bytes y MD5, seguida de
  validación completa de gzip y estructura FASTQ.
- Entorno bioinformático `cacao-microbiome` en Ubuntu/WSL2 y verificador que
  ejecuta las herramientas de consola, carga los paquetes R y registra 18
  versiones resueltas en una tabla trazable y un lock explícito `linux-64`.
- Target Snakemake de QC crudo con 12 trabajos FastQC, agregación MultiQC y tabla
  validada que cruza métricas con la integridad y los conteos FASTQ previos.
- Configuración exacta y trazable de los primers 515F/806R y sus constructos
  L1/L2 para `PRJNA492720`, verificados en el suplemento primario Table S2.
- Target Snakemake y tabla machine-readable con 96 pruebas Cutadapt de presencia,
  orientación y constructo completo sobre los 12 FASTQ del piloto.
- Recorte pareado y configurable de primers del piloto con Cutadapt, sin
  descarte de lecturas sin primer, y resumen validado por dirección.

### Changed

- Nombre y descripción inicial del proyecto actualizados a
  **Cacao Fermentation Microbiome**.
- `PRJNA492720` adoptado como piloto con 94 corridas 16S V4 compatibles con la
  publicación; dos estudios permanecen pendientes y tres WGS se excluyen.
- `PRJNA627078` incorporado como segundo estudio primario con 60 corridas 16S
  V3-V4 longitudinales; tres proyectos prioritarios quedan pendientes y
  `PRJNA962540` se excluye por ausencia de una serie temporal 16S.
- `PRJNA492720` confirmado objetivamente como piloto; su vertical slice de 12
  FASTQ pasó checksum e integridad estructural sin versionar los datos crudos.
- Canales del entorno restringidos a conda-forge/Bioconda con prioridad estricta
  y restricciones de ABI explícitas para DADA2, ANCOMBC, Matrix y lme4.
