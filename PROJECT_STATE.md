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
- Seis estudios cribados en detalle: uno incluido, dos pendientes y tres
  excluidos de Fase I.
- Inventario generado: 861 corridas, 131 muestras candidatas, 735 exclusiones
  trazables y 94 corridas incluidas en el piloto `PRJNA492720`.
- Validador independiente y 14 tests unitarios aprobados.

## En progreso

- Revisión de publicación y elegibilidad de 18 BioProjects detectados por la
  consulta sistemática; cinco tienen prioridad por señales simultáneas de
  amplicón, 16S/bacterias y temporalidad.
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

1. Cribar los cinco proyectos de revisión prioritaria y localizar sus artículos.
2. Resolver los otros 13 candidatos de revisión manual y registrar cada decisión.
3. Congelar el inventario elegible de Fase I y generar el manifest de descarga.
4. Ejecutar una descarga piloto contenida de `PRJNA492720` con checksums.
5. Derivar umbrales de QC desde el piloto antes de automatizar DADA2.

## Limitaciones conocidas

- No se instaló el entorno bioinformático completo: en este equipo están
  ausentes R, Snakemake, FastQC, MultiQC y Cutadapt.
- No se descargaron FASTQ ni se generaron resultados biológicos.
- La revisión paper-level de 18 candidatos todavía no está cerrada.
- Los metadatos de primers y diseño temporal son incompletos en varios depósitos.
