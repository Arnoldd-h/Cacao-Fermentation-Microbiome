# Estado del proyecto

Última actualización: 2026-10-06 (ejecuciones del 2026-10-07 UTC).

## Fase actual

Inventario de 34 candidatos consolidado y QC reproducido con procedencia.
Piloto DADA2 implementado y primera ejecución técnica validada; falta cerrar
su ejecución integrada final de Snakemake con el código versionado.

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
- Primera ejecución DADA2: 165.034 pares filtrados, 160.649 ensamblados y 156.561
  sin quimeras (88,55 % de la entrada), 93 ASVs y seis muestras conservadas.
- Ambos modelos de error convergieron en cinco rondas; gráficas revisadas y
  doce warnings de representación de frecuencias cero preservados.
- Validador independiente de DADA2: conteos entre etapas, tabla ASV/FASTA,
  distribuciones de longitud y hashes; regresiones Python y pruebas R.

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

## En progreso

- Ejecución final integrada DADA2 después de congelar el código y el QC.
- Publicación local de los hitos Git y cierre de documentación.

El inventario pasó el validador independiente; las unidades analíticas cubren
182 runs, 95 observaciones lote-tiempo y 12 lotes. QC terminó 45 trabajos y
los 24 FASTQ raw/interim coinciden byte a byte con los archivos previos.

## Bloqueos y limitaciones

- Siete estudios permanecen pendientes por falta o inconsistencia de evidencia
  pública; no se inventan primers, tiempos, tratamientos ni asignaciones de run.
- El procesamiento real abarca seis runs del piloto; el resto de los FASTQ no
  se ha descargado. No hay aún taxonomía, diversidad ni inferencia biológica.
- Las ASVs de 237–243 nt no se filtran por una longitud canónica; el depósito
  tiene evidencia de recorte previo y su procesamiento exacto es incompleto.
- Los residuos raros en orientación inversa y la ausencia de controles negativos
  requieren interpretación posterior; no demuestran contaminación ni su ausencia.
- País, estudio, tecnología y región 16S pueden estar confundidos; no se les
  atribuirán efectos separables sin respaldo del diseño.
- Licencia y metadatos de autoría del repositorio siguen pendientes del titular.

## Próximas tareas

1. Cerrar validaciones y tags de inventario, QC y piloto DADA2.
2. Fijar versión y checksum de una referencia taxonómica y validar la clasificación.
3. Escalar el estudio piloto completo conservando parámetros por estudio y
   revisar por separado los pipelines V3-V4 y PacBio antes de integrar taxonomía.
4. Definir tratamiento de ceros, agregación de submuestras y modelos longitudinales
   antes de estimar efectos, con sensibilidad temporal absoluta y relativa.

## Control de versiones

Trabajo en `codex/pilot-dada2-reproducibility`; `main` permanece disponible.
Los commits son locales. No se publica ni se hace push sin autorización.
Los datos pesados permanecen ignorados; las salidas pequeñas conservan procedencia.
El commit `44f5f4a` prerregistra los parámetros; `eafadd7` corrige identificadores
SRA; `a2bfa82` y `dd4f364` cierran reproducibilidad del workflow QC.
Para el checkout exacto usar `git log -1 --oneline` y `git status`.
