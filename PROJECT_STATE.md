# Estado del proyecto

Última actualización: 2026-10-09.

## Fase actual

Inventario, QC y piloto DADA2 consolidados e integrados en `main` mediante PR #1
(`9f0a478`). Clasificación taxonómica del piloto ejecutada y validada con
SILVA 138.2, checksum de origen y umbrales fijados antes de observar las
asignaciones. Se conservan 93 ASVs y 156.561 lecturas de seis muestras.
Taxonomía integrada en `main` mediante PR #2 (`6da8f86`) y tag `v0.4.0-taxonomy`.
Separación auditable y diversidad descriptiva del piloto ejecutadas y validadas:
88 ASVs y 123.693 lecturas derivadas, con las seis muestras conservadas.
Los criterios, métricas y sustitución de ceros quedaron registrados antes de
calcular diversidad. El hito `v0.5.0-diversity` se limita a esta descripción;
la inferencia temporal y la integración entre estudios siguen pendientes.
La diversidad está integrada en `main` mediante PR #3 (`868afa8`). En progreso:
ampliación al estudio completo de 94 runs, con auditoría de fuentes, recursos
limitados y conservación de los parámetros científicos del piloto.
NCBI/ENA confirman los 94 runs incluidos. La auditoría de FASTQ identifica
93 runs con archivos pareados y uno sin pares completos (SRR7899803); su exclusión de la
ruta paired-end se registra aparte y no cambia la inclusión del inventario.

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
- Suite de la fase DADA2: 91 pruebas Python en Linux aprobadas, incluidas seis integraciones
  Snakemake; 32 comprobaciones R y regresión pareada SRA aprobadas.
- Corrida final del 2026-10-09 UTC (8 de octubre en Colombia), sobre `6221445`:
  18 entradas y 18 artefactos DADA2 validados. Las cuatro tablas principales
  son idénticas a las de la primera corrida técnica.
- SILVA 138.2 descargada y verificada por bytes, MD5, SHA-256 y lectura completa
  de gzip/FASTA: 452.055 secuencias; base protegida y fuera de Git.
- Primera clasificación del 2026-10-09 sobre `d13c9eb`, con checkout limpio al iniciar:
  21 entradas y 15 artefactos pasan validación independiente, sin warnings R.
- Bootstrap 80: 64/93 ASVs asignadas a género, asociadas a 99.481 lecturas
  (63,54 %). Sensibilidad 50: 79/93 ASVs y 104.974 lecturas (67,05 %).
- Cinco ASVs marcadas como Chloroplast/Mitochondria, con 32.868 lecturas
  (20,99 %), conservadas sin exclusión automática. Tablas y figuras trazables.
- Suite de la fase taxonómica: 111 pruebas Python en Linux y seis comprobaciones taxonómicas
  R aprobadas; dry-run integrado `all` sin trabajos pendientes.
- Repetición sobre `4e7ce24` tras ajustar la leyenda: siete tablas científicas
  idénticas por SHA-256, 21 entradas y 15 artefactos validados. Figuras revisadas.
- Política de separación y diversidad registrada en `f9132b5`; implementación
  sobre `5221b0f`, conservando las tablas originales y sin nuevas dependencias.
- Taxonomía renovada sobre `5221b0f`: siete tablas científicas idénticas por
  SHA-256 a las del hito anterior, con validación de 21 entradas/15 artefactos.
- Separación de cinco ASVs de orgánulos (32.868 lecturas), con registro de las
  93 ASVs y balance por muestra: 88 ASVs y 123.693 lecturas retenidas.
  Validador independiente de 21 entradas/nueve artefactos; ninguna muestra excluida.
- Diversidad alfa sin rarefacción, CLR/Aitchison con pseudoconteos 1 y 0,5,
  Bray–Curtis sobre proporciones y PCA: 20 entradas/17 artefactos validados.
  Profundidad 9.250–36.508 lecturas; riqueza observada 2–46 ASVs.
- Figuras alfa/PCA en PDF/SVG/PNG de 300 dpi, manteniendo muestras individuales
  y escala geométrica igual en los ejes PCA; sin warnings R de diversidad.
- Suite actual: 125 pruebas Python en Linux aprobadas. Las 12 específicas
  de diversidad también pasan tras los ajustes de figuras.
- Corrida final de diversidad sobre `342dd4e`: seis tablas científicas y
  versiones idénticas por SHA-256 a la primera ejecución, figuras revisadas
  y 20 entradas/17 artefactos validados. Dry-run integrado sin trabajos pendientes.

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
Separación y diversidad tienen productores y validadores propios, con política
en `config/diversity.yaml`. Las tablas derivadas están en `data/processed/`.
Cambiar el pseudoconteo invalida esas etapas sin reclasificar taxonomía;
esta conducta se comprueba con una regresión del workflow.
El dry-run informa que 12 descargas previas carecen de metadata histórica del
scheduler Snakemake. Sus FASTQ siguen verificados por manifest, bytes, MD5,
validación estructural y procedencia científica; no se modifican los crudos
para reconstruir ese historial auxiliar.
No queda una ejecución activa.

## Bloqueos y limitaciones

- Siete estudios permanecen pendientes por falta o inconsistencia de evidencia
  pública; no se inventan primers, tiempos, tratamientos ni asignaciones de run.
- El procesamiento real abarca seis runs del piloto; el resto de los FASTQ no
  se ha descargado. Hay taxonomía y diversidad descriptiva del piloto; la
  inferencia biológica sigue pendiente. Seis runs de cuatro lotes y distintos
  estratos no reconstruyen sus trayectorias ni equivalen a seis lotes independientes.
- El 20,99 % de las lecturas pertenece a ASVs con etiquetas de orgánulos.
  La cobertura taxonómica actual incluye esas lecturas en su denominador.
  Su separación está configurada y registrada en la tabla derivada; las
  etiquetas no identifican especie huésped ni demuestran contaminación.
- Profundidad desigual y dependencia del pseudoconteo respecto a la escala de
  conteos limitan las comparaciones descriptivas. La riqueza no está rarefactada;
  no se calcula significación ni se elige la sensibilidad por separación visual.
- Las ASVs de 237–243 nt no se filtran por una longitud canónica; el depósito
  tiene evidencia de recorte previo y su procesamiento exacto es incompleto.
- Los residuos raros en orientación inversa y la ausencia de controles negativos
  requieren interpretación posterior; no demuestran contaminación ni su ausencia.
- País, estudio, tecnología y región 16S pueden estar confundidos; no se les
  atribuirán efectos separables sin respaldo del diseño.
- Licencia y metadatos de autoría del repositorio siguen pendientes del titular.

## Próximas tareas

1. Escalar los 94 runs del estudio piloto, con estimación de recursos,
   reconstrucción de FASTQ y validación por estudio antes de interpretar efectos.
2. Conservar parámetros por estudio y
   revisar por separado los pipelines V3-V4 y PacBio antes de integrar taxonomía.
3. Definir tratamiento de ceros, agregación de submuestras y modelos longitudinales
   antes de estimar efectos, con sensibilidad temporal absoluta y relativa.

## Control de versiones

La etapa actual se desarrolló en `codex/pilot-diversity`, creada desde `main`
con taxonomía integrada (`6da8f86`). Su hito es `v0.5.0-diversity`, limitado
a separación y diversidad descriptiva del piloto. `f9132b5` registra la
política; `5221b0f` implementa el workflow y ejecuta taxonomía/separación;
`dc1cf03` ajusta las figuras. El commit exacto de la ejecución final de
diversidad queda en su procedencia y marcador `SUCCESS`.
Los datos pesados permanecen ignorados; las salidas pequeñas conservan procedencia.
El commit `44f5f4a` prerregistra los parámetros; `eafadd7` corrige identificadores
SRA; `a2bfa82` y `dd4f364` cierran reproducibilidad del workflow QC.
`20feeea` corrige reevaluación del checkpoint; `32ab045` consolida inventario y
`6221445` consolida QC. Los hitos publicados son `v0.1.0-dataset-inventory`,
`v0.2.0-pilot-qc`, `v0.3.0-pilot-dada2` y `v0.4.0-taxonomy`.
La procedencia DADA2 registra `git_dirty=true` por documentación y resultados
previos todavía sin commit al iniciar. Sus seis entradas versionadas (config,
manifest, scripts R y tablas QC) coinciden byte a byte con `6221445`.
La primera taxonomía registra `git_dirty=false` sobre `d13c9eb`. Las salidas
del hito taxonómico proceden de `4e7ce24`: antes de lanzar el workflow el checkout estaba
limpio; Snakemake retiró sus salidas previas para regenerarlas y el ejecutor
registra `git_dirty=true` sólo por esas eliminaciones. El código y las entradas
científicas permanecen fijados por sus hashes; las siete tablas reproducen la
primera corrida sin cambios.
La renovación taxonómica actual procede de `5221b0f` y conserva esas siete
tablas. Separación y diversidad registran también los resultados anteriores
o recién generados que todavía no tenían commit al iniciar. Ese estado Git
se conserva íntegro, junto a los hashes de todas las entradas y salidas.
Para el checkout exacto usar `git log -1 --oneline` y `git status`.
