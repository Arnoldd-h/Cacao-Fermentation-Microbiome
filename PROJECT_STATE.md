# Estado del proyecto

Última actualización: 2026-10-09.

## Fase actual

Inventario, QC y piloto DADA2 consolidados e integrados en `main` mediante PR #1
(`9f0a478`). Clasificación taxonómica del piloto ejecutada y validada con
SILVA 138.2, checksum de origen y umbrales fijados antes de observar las
asignaciones. Se conservan 93 ASVs y 156.561 lecturas de seis muestras.
La siguiente etapa requiere separar de forma auditable las etiquetas de
orgánulos antes de preparar diversidad y ampliar el estudio completo.

## Completado

- Configuración, protocolo, registro de decisiones, inventario y validadores.
- Revisión de los 34 BioProjects descubiertos: cuatro incluidos, siete pendientes
  por evidencia insuficiente y 23 excluidos; 2.357 runs revisados y 182 incluidos.
- Selección del piloto PRJNA492720 y manifest de seis runs, sin cambios tras
  ampliar el inventario; 12 FASTQ crudos con bytes, MD5, gzip y estructura válidos.
- Entorno Linux/WSL2 verificado: FastQC 0.12.1, MultiQC 1.25.1, Cutadapt 5.0,
  Snakemake 8.30.0, R 4.4.3 y DADA2 1.34.0; 18 componentes y lock explícito.
- Ejecutor Micromamba que resuelve el prefijo registrado, sin rutas absolutas
  codificadas ni dependencia de una raíz Micromamba particular.
- Workflow reconstruible sin manifest/FASTQ preexistentes, protección de datos
  crudos y comparación por contenido del manifest para evitar repetir QC.
- Validación de parámetros realmente ejecutados por Cutadapt y procedencia de
  resúmenes QC con hashes de entradas, código, configuración y salidas.
- Jerarquía run/muestra/lote/tiempo: cada lote es una unidad independiente;
  los tiempos son medidas repetidas y las submuestras permanecen identificadas.
- QC inicial completo: 176.798 pares, 100 % de retención Cutadapt, siete R2
  recortados y 238 bases retiradas; calidad por base y adaptadores pasan 12/12.
- DADA2 configurado antes de ejecutarse; corrección del parser SRA manteniendo
  la comprobación obligatoria de identificadores pareados.
- Ejecución integrada DADA2: 165.034 pares filtrados, 160.649 ensamblados y 156.561
  sin quimeras (88,55 % de la entrada), 93 ASVs y seis muestras conservadas.
- Ambos modelos de error convergieron en cinco rondas; gráficas revisadas y
  doce warnings de representación de frecuencias cero preservados.
- Validador independiente de DADA2: conteos entre etapas, tabla ASV/FASTA,
  distribuciones de longitud y hashes; regresiones Python y pruebas R.
- Suite final: 91 pruebas Python en Linux aprobadas, incluidas seis integraciones
  Snakemake; 32 comprobaciones R y regresión pareada SRA aprobadas.
- Corrida final del 2026-10-09 UTC (8 de octubre en Colombia), sobre `6221445`:
  18 entradas y 18 artefactos DADA2 validados. Las cuatro tablas principales
  son idénticas a las de la primera corrida técnica.
- SILVA 138.2 descargada y verificada por bytes, MD5, SHA-256 y lectura completa
  de gzip/FASTA: 452.055 secuencias; base protegida y fuera de Git.
- Clasificación del 2026-10-09 sobre `d13c9eb`, con checkout limpio al iniciar:
  21 entradas y 15 artefactos pasan validación independiente, sin warnings R.
- Bootstrap 80: 64/93 ASVs asignadas a género, asociadas a 99.481 lecturas
  (63,54 %). Sensibilidad 50: 79/93 ASVs y 104.974 lecturas (67,05 %).
- Cinco ASVs marcadas como Chloroplast/Mitochondria, con 32.868 lecturas
  (20,99 %), conservadas sin exclusión automática. Tablas y figuras trazables.
- Suite actual: 111 pruebas Python en Linux y seis comprobaciones taxonómicas
  R aprobadas; dry-run integrado `all` sin trabajos pendientes.

## Datasets incluidos

| BioProject | Runs incluidos | Lotes | Tecnología |
|---|---:|---:|---|
| PRJNA492720 | 94 | 4 | Illumina paired-end, V4; piloto |
| PRJNA627078 | 60 | 4 | Illumina paired-end, V3-V4 |
| PRJEB40850 | 12 | 2 | Illumina paired-end, V4 |
| PRJEB57747 | 16 | 2 | PacBio Sequel II, full-length, single-end |

Se incluyen 182 runs y 12 lotes, sin confundir runs con réplicas independientes.
PacBio requiere un procesamiento específico y no entra en el piloto Illumina.
Las fuentes y decisiones por estudio están en `config/datasets.yaml` y
`reports/inventory_review_2026-10-04.md`.

## Estado operativo

El inventario pasó el validador independiente; las unidades analíticas cubren
182 runs, 95 observaciones lote-tiempo y 12 lotes. QC terminó 45 trabajos y
los 24 FASTQ raw/interim coinciden byte a byte con los archivos previos.
La corrección de los consumidores del checkpoint elimina la repetición de
resúmenes: el dry-run conjunto QC termina sin trabajos pendientes.
Los cinco registros de procedencia renovados pasan 222 comprobaciones de
bytes y SHA-256 de configuración, código, entradas y salidas.
El target `all` completó descubrimiento, DADA2 y validación; el siguiente
dry-run devuelve `Nothing to be done`. El descubrimiento actualizado conserva
34 BioProjects y sus conteos de runs coinciden con el inventario versionado.
La taxonomía añade validación independiente de soportes, máscaras jerárquicas,
marcas y cobertura por muestra; preserva las tablas y configuración DADA2.
No queda una ejecución activa.

## Bloqueos y limitaciones

- Siete estudios permanecen pendientes por falta o inconsistencia de evidencia
  pública; no se inventan primers, tiempos, tratamientos ni asignaciones de run.
- El procesamiento real abarca seis runs del piloto; el resto de los FASTQ no
  se ha descargado. Hay taxonomía del piloto; diversidad e inferencia biológica
  siguen pendientes.
- El 20,99 % de las lecturas pertenece a ASVs con etiquetas de orgánulos.
  La cobertura taxonómica actual incluye esas lecturas en su denominador.
  Su separación posterior exige configuración, decisión y registro explícitos;
  estas etiquetas no identifican especie huésped ni demuestran contaminación.
- Las ASVs de 237–243 nt no se filtran por una longitud canónica; el depósito
  tiene evidencia de recorte previo y su procesamiento exacto es incompleto.
- Los residuos raros en orientación inversa y la ausencia de controles negativos
  requieren interpretación posterior; no demuestran contaminación ni su ausencia.
- País, estudio, tecnología y región 16S pueden estar confundidos; no se les
  atribuirán efectos separables sin respaldo del diseño.
- Licencia y metadatos de autoría del repositorio siguen pendientes del titular.

## Próximas tareas

1. Registrar la política de separación de orgánulos y generar conteos derivados
   con un registro por ASV, conservando las tablas originales y sin elegir
   exclusiones según efectos temporales.
2. Escalar el estudio piloto completo conservando parámetros por estudio y
   revisar por separado los pipelines V3-V4 y PacBio antes de integrar taxonomía.
3. Definir tratamiento de ceros, agregación de submuestras y modelos longitudinales
   antes de estimar efectos, con sensibilidad temporal absoluta y relativa.

## Control de versiones

Trabajo en `codex/pilot-taxonomy`, creada desde el `main` integrado (`9f0a478`).
La fase anterior y sus tres tags están publicados. Los nuevos commits son
locales hasta publicar la siguiente unidad completa dentro de la autorización.
Los datos pesados permanecen ignorados; las salidas pequeñas conservan procedencia.
El commit `44f5f4a` prerregistra los parámetros; `eafadd7` corrige identificadores
SRA; `a2bfa82` y `dd4f364` cierran reproducibilidad del workflow QC.
`20feeea` corrige reevaluación del checkpoint; `32ab045` consolida inventario y
`6221445` consolida QC. Los hitos publicados son `v0.1.0-dataset-inventory`,
`v0.2.0-pilot-qc` y `v0.3.0-pilot-dada2`.
La procedencia DADA2 registra `git_dirty=true` por documentación y resultados
previos todavía sin commit al iniciar. Sus seis entradas versionadas (config,
manifest, scripts R y tablas QC) coinciden byte a byte con `6221445`.
La procedencia taxonómica registra `git_dirty=false` y el commit `d13c9eb`.
Para el checkout exacto usar `git log -1 --oneline` y `git status`.
