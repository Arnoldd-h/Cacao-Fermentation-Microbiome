# Estado del proyecto

Última actualización: 2026-08-19

## Fase actual

Inventario sistemático de datasets públicos y cribado científico de candidatos
para la Fase I bacteriana 16S.

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
- Validador independiente y 17 tests unitarios aprobados.

## En progreso

- Revisión de publicación y elegibilidad de 13 BioProjects detectados por la
  consulta sistemática que todavía no están configurados.
- Verificación de primers, región 16S, diseño temporal y condición espontánea
  donde la metadata pública es insuficiente.

## Bloqueado

- Sin bloqueos registrados.

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

## Próximas tareas

1. Generar y validar la comparación reproducible del dataset piloto.
2. Resolver los 13 candidatos de revisión manual y registrar cada decisión.
3. Congelar el inventario elegible de Fase I y generar el manifest de descarga.
4. Ejecutar una descarga piloto contenida de `PRJNA492720` con checksums.
5. Derivar umbrales de QC desde el piloto antes de automatizar DADA2.

## Limitaciones conocidas

- No se instaló el entorno bioinformático completo: en este equipo están
  ausentes R, Snakemake, FastQC, MultiQC y Cutadapt.
- No se descargaron FASTQ ni se generaron resultados biológicos.
- La revisión paper-level de 13 candidatos todavía no está cerrada.
- Los metadatos de primers y diseño temporal son incompletos en varios depósitos.
