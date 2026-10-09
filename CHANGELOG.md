# Changelog

Los cambios importantes de este proyecto se documentarán en este archivo.

## [Unreleased]

### Added

- Pipeline DADA2 por estudio para el piloto Illumina: filtros pareados,
  aprendizaje de errores, inferencia, ensamblaje, quimeras y diagnósticos.
- Validación independiente de conservación de lecturas, ASVs, FASTA,
  distribución de longitudes y checksums de entradas y salidas de DADA2.
- Pruebas R de contratos y regresión de identificadores SRA; target
  `pilot_dada2` integrado en el workflow principal.
- Workflow del piloto reconstruible mediante checkpoint del manifest y reglas
  productoras de FASTQ protegidos; pruebas de arranque sin manifest ni datos.
- Validación de parámetros y rutas realmente ejecutados por Cutadapt, selección
  de primers por BioProject y procedencia de todos los resúmenes QC.
- Jerarquía reproducible run/muestra/lote/tiempo que conserva las submuestras
  anidadas y evita contarlas como fermentaciones independientes.
- Registros de procedencia con fecha, commit, estado del checkout, parámetros
  y checksums de entradas, configuración y salidas.
- Ejecutor portable de Micromamba que descubre el prefijo registrado o acepta
  una selección local explícita y rechaza instalaciones ausentes o ambiguas.
- Parámetros iniciales de DADA2 por estudio fijados antes de la ejecución del
  piloto, con semilla, filtros, aprendizaje, merging y quimeras explícitos.
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
- QC post-Cutadapt reproducible con 12 FastQC, MultiQC y comparación tabular
  contra el QC crudo; se registran retención, longitudes, estados FastQC y
  procedencia de cada dirección.
- Detección de primers post-Cutadapt integrada en el workflow y objetivo final
  de Snakemake actualizado para exigir todos los artefactos del QC del piloto.

### Changed

- Los consumidores y agregados QC/DADA2 resuelven directamente el checkpoint
  del manifest; desaparecen razones de actualización obsoletas que repetían
  resúmenes pese a descartar sus productores. Se prueba estabilidad del grafo.
- Revisión de los 34 BioProjects: cuatro incluidos (182 runs), siete pendientes
  (93 runs candidatos) y 23 excluidos. Se incorporan los controles espontáneos
  de Costa Rica 2017/2019, conservando PacBio separado del pipeline Illumina.
- PRJNA865318 y PRJNA1104253 pasan de pendientes a excluidos por evidencia
  primaria; cada pendiente restante tiene un requisito de resolución explícito.
- Metadata de archivos enviados permite separar 16S/ITS sin descargar lecturas;
  lotes, duración y campos no verificados permanecen vacíos en vez de inferirse.
- La reescritura idéntica de tablas conserva su archivo y el workflow compara
  el contenido SHA-256 del manifest, evitando repetir QC por un cambio de fecha.
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
