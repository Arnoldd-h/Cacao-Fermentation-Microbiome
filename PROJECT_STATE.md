# Estado del proyecto

Última actualización: 2026-08-19

## Fase actual

Inventario sistemático de datasets públicos para la Fase I bacteriana 16S y
vertical slice técnica del piloto.

## Completado

- Repositorio remoto renombrado a `Cacao-Fermentation-Microbiome`.
- Rama estable `main` clonada y conectada con `origin`.
- Política de commits, trazabilidad y manejo de datos incorporada al repositorio.
- Exclusiones iniciales para datos ómicos pesados, credenciales, entornos y cachés.
- Pregunta, hipótesis, protocolo, criterios de inclusión/exclusión y decisiones
  metodológicas iniciales documentados.
- Arquitectura reproducible con configuración externa, librería Python,
  scripts, workflow Snakemake, entorno declarado y tests.
- Skill local `sra-metadata-inventory` creado y validado para regenerar y auditar
  la metadata sin descargar FASTQ.
- Consulta sistemática NCBI/ENA ejecutada: 34 BioProjects descubiertos.
- Once estudios cribados en detalle: dos incluidos, cinco pendientes y cuatro
  excluidos de Fase I.
- Inventario generado: 1.573 corridas, 255 muestras candidatas, 1.327
  exclusiones trazables y 154 corridas incluidas.
- `PRJNA627078` incorporado como segundo estudio primario con 60 corridas V3-V4
  tradicionales y cinco tiempos.
- Comparación reproducible de 11 estudios: `PRJNA492720` seleccionado como
  piloto frente a `PRJNA627078`.
- Manifest equilibrado de seis runs (dos por etapa temporal y lotes distintos)
  con 25.154.687 bytes estimados.
- Doce FASTQ descargados localmente: 12/12 bytes, MD5, gzip y estructura FASTQ
  válidos; 353.596 reads y 78.651.325 bases en total.
- Validador independiente y 30 tests unitarios aprobados.

## Datasets incluidos

- `PRJNA492720`: 94 runs 16S V4; piloto seleccionado y seis runs descargados.
- `PRJNA627078`: 60 runs 16S V3-V4 tradicionales.

## Datasets pendientes

- `PRJNA865318`, `PRJNA1104253`, `PRJEB53853`, `PRJEB82327` y `PRJEB82871`.

## Datasets excluidos de Fase I

- `PRJNA552479`, `PRJNA1257864`, `PRJNA1264670` y `PRJNA962540`.

## En progreso

- Revisión de publicación y elegibilidad de 13 BioProjects detectados por la
  consulta sistemática que todavía no están configurados.
- Verificación de primers, región 16S, diseño temporal y condición espontánea
  donde la metadata pública es insuficiente.
- Preparación del entorno bioinformático para FastQC, MultiQC y Cutadapt.

## Bloqueado

- El QC con FastQC/MultiQC y el procesamiento posterior requieren instalar el
  entorno declarado; no hay Conda, micromamba, R ni las herramientas del pipeline
  disponibles en esta máquina.

## Decisiones importantes

- `main` debe conservar un estado razonablemente funcional y reproducible.
- Los datos crudos no se versionarán; se reconstruirán desde accessions,
  manifests y scripts.
- Los cambios metodológicos deberán quedar documentados y validados.
- Los commits locales se crearán por unidades lógicas; los pushes requerirán
  autorización explícita.
- La triaje automática no puede promover proyectos a inclusión: toda inclusión
  exige revisión científica documentada en `config/datasets.yaml`.
- `PRJNA492720` es el piloto inicial; sus 94 corridas paper-consistent se
  seleccionan por reglas explícitas y auditables.
- La vertical slice limita la descarga a dos muestras por etapa y valida cada
  archivo antes de cualquier QC o recorte.

## Próximas tareas

1. Resolver los 13 candidatos de revisión manual y congelar el inventario de Fase I.
2. Instalar y verificar el entorno bioinformático reproducible.
3. Ejecutar FastQC/MultiQC y detección de primers sobre los seis runs.
4. Definir Cutadapt desde primers verificados y perfiles de calidad reales.
5. Probar DADA2 por estudio y derivar parámetros antes de escalar a 94 runs.

## Limitaciones conocidas

- No se instaló el entorno bioinformático completo: están ausentes R, Snakemake,
  FastQC, MultiQC y Cutadapt.
- Solo se descargó la vertical slice de 12 FASTQ; no se descargó el estudio completo.
- No se ejecutaron FastQC, recorte, DADA2, taxonomía ni resultados biológicos.
- La revisión paper-level de 13 candidatos todavía no está cerrada.
- Los metadatos de primers y diseño temporal son incompletos en varios depósitos.

## Último commit relevante

`c1ecf57 feat: add resumable pilot FASTQ validation`, usado para generar los
reportes QC de descarga e integridad.
