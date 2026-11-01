# Estado del proyecto

Última actualización: 2026-10-31.

## Fase actual

Control explícito del desempate taxonómico validado en el entorno registrado:
adaptador `r_seeded_query_mt19937_v1`, sin modificar el paquete DADA2 instalado.
Piloto y estudio completo repiten exactamente sus siete tablas científicas
entre procesos independientes con dos y un hilo. Método registrado antes
de reclasificar en `e008f2b`; tablas históricas completas archivadas.

Primer estudio completo procesado y validado de forma descriptiva: PRJNA492720,
93 muestras pareadas de cuatro lotes y 47 observaciones lote-tiempo. El
inventario conserva 94 runs; SRR7899803 no ofrece pares completos y se excluye
solo de esta ruta, con decisión auditable. El bloqueo de ENA quedó resuelto:
186/186 FASTQ, 394.515.770 bytes, con bytes, MD5, gzip/CRC y estructura válidos.

QC conserva 2.728.920 pares tras Cutadapt y pasa 505 entradas/11 artefactos.
DADA2 conserva 782 ASVs y 2.254.530 pares sin quimeras (82,62 %), con
192 entradas/18 artefactos comprobados. Ambos modelos convergen en seis rondas.
SILVA 138.2 clasifica género en 475 ASVs a bootstrap 80; se validan 26
entradas/15 artefactos. La separación de 18 ASVs etiquetadas como orgánulos
conserva 764 ASVs, 2.096.980 lecturas y todas las muestras. Separación y
diversidad pasan 21 entradas/nueve artefactos y 20 entradas/17 artefactos.

Los resúmenes por familia/género contienen 1.477 filas, sin imputar etiquetas
ni renormalizar solo taxones clasificados (15 entradas/cinco artefactos).
Las figuras temporales representan 93 muestras y ocho series de lote/estrato,
con horas originales y sin conectar huecos; exportación PDF/SVG/PNG a 300 dpi,
16 entradas/ocho artefactos y revisión visual. Las fuentes exactas quedan en
la procedencia de cada etapa y en el registro de ejecución del 2026-10-31.

El primer recorrido cerró 481 trabajos; la renovación de la figura taxonómica
cerró siete trabajos adicionales sin repetir QC ni DADA2 completos. La suite
completa pasa 161 pruebas Python en Linux (135,991 s). El piloto de seis muestras
cerró 33 trabajos con el código compartido actual: ocho tablas DADA2, ocho
taxonómicas y once derivadas/de diversidad coinciden con la referencia anterior.
Sus cinco registros QC pasan 257 comprobaciones de bytes/SHA-256. Es un
subconjunto del estudio completo, no un
segundo estudio ni seis muestras independientes adicionales.

La comparación del adaptador con el histórico conserva taxonomía primaria,
44 tablas TSV de DADA2/separación/diversidad y ocho artefactos TSV/PNG de
composición/tiempo. Sólo ASV000100 (cinco lecturas, BactSaI012z1) cambia su
soporte de familia/género de 48 a 51: entra en la sensibilidad a 50 como
Yersiniaceae/Samsonia, pero permanece sin esas etiquetas a 80. Género a 50:
649 ASVs/2.051.014 lecturas. El piloto conserva sus siete tablas históricas.
No se cambia un umbral ni se selecciona la corrida por sus resultados.
El contrato exige los mismos datos/orden y entorno; no garantiza igualdad
entre plataformas, compiladores o al reordenar entradas.

Validación de esta continuación: 166 pruebas de suite (140,517 s), más una
regresión adicional de auditoría que rechaza entradas discordantes aunque
las tablas coincidan. Los siete trabajos de cada alcance y las repeticiones
independientes pasan sus validadores; la cobertura renovada se revisó visualmente.

La inferencia longitudinal, el core temporal y la integración/meta-análisis
siguen pendientes. Los otros tres estudios incluidos todavía no están procesados.
La ampliación se entrega mediante PR #6; PR #5 integró recuperación de datos,
recursos y resúmenes de composición en `main` (`7cd3986`). Los FASTQ, RDS,
base SILVA e informes pesados permanecen fuera de Git.

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
- Suite del hito piloto: 125 pruebas Python en Linux aprobadas. Las 12 específicas
  de diversidad también pasan tras los ajustes de figuras.
- Corrida final de diversidad sobre `342dd4e`: seis tablas científicas y
  versiones idénticas por SHA-256 a la primera ejecución, figuras revisadas
  y 20 entradas/17 artefactos validados. Dry-run integrado sin trabajos pendientes.
- Workflow de estudio completo sobre `0cc8140`, con corrección de referencia
  compartida `fa236ed`, decisiones de elegibilidad `b10a15e` y alcance `8fd5334`.
- Manifest de 93 runs/186 archivos a partir de 94 runs auditados; los 189 streams
  están registrados y tres streams sin pareja quedan excluidos explícitamente.
- Auditorías de manifest/fuentes/empaquetado: 31 entradas y seis artefactos del
  manifest, ocho registros de fuentes y 14 de empaquetado pasan hashes actuales.
- Suite final: 144 pruebas Python en Linux aprobadas, incluidas 19 de ampliación.
- Piloto renovado sobre `e472984`: siete tablas TSV idénticas por SHA-256,
  20 entradas/17 artefactos validados y dry-run estable sin trabajos pendientes.
- Recuperación del R2 pendiente desde ENA el 2026-10-31 y cierre de descarga:
  186/186 FASTQ, 394.515.770 bytes; cero exclusiones añadidas por transporte.
- Entorno existente verificado de nuevo: 18 componentes, con las mismas versiones;
  cuatro regresiones de descarga aprobadas y dry-run de ampliación de 670 trabajos.
- Suite completa renovada: 144 pruebas aprobadas. Ajuste operativo a ocho tareas
  concurrentes y 4.500 MB, conservando dos hilos para DADA2/taxonomía y los métodos.
- Resúmenes de composición del piloto: 95 filas por familia/género sobre seis
  muestras, 88 ASVs y 123.693 lecturas; 15 entradas/cinco artefactos validados,
  sin renormalizar sólo taxones clasificados. Siete pruebas nuevas aprobadas.
- Suite completa posterior a la implementación de composición: 151 pruebas
  Python aprobadas en Linux, incluidas las integraciones R y Snakemake existentes.
- Corrección de rutas Cutadapt en `67c3820`: 154 pruebas completas aprobadas,
  incluidos casos de rutas completas, escapes de directorio y una integración real
  con FASTQ sintéticos. Las notas del manuscrito documentan la evidencia preliminar
  y las hipótesis pendientes sin presentar el análisis completo como terminado.
- QC completo de 93 runs: 2.728.920 pares, retención Cutadapt 100 %, tres R1 y
  109 R2 recortados, 3.812 bases retiradas; 505 entradas/11 artefactos validados.
  Calidad por base: 184/186 pasan antes y después del recorte; se conservan los
  avisos de una muestra y las alertas de adaptadores de dos muestras, sin exclusiones.

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

## Estado operativo y antecedentes

En el hito piloto anterior, el inventario pasó el validador independiente; las unidades analíticas cubren
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
La descarga completa está validada desde el 2026-10-31. El archivo HTML de
797 bytes del intento anterior se conserva en `results/tmp/`, fuera de Git;
no se utilizó para reanudar bytes de un FASTQ. QC, DADA2, taxonomía, separación,
diversidad, composición y figuras temporales de los 93 runs están validados.
El piloto cerró 33 trabajos tras cambios compartidos de QC y del título
taxonómico, manteniendo sus tablas. La normalización local a LF del config
completo renovó cinco trabajos de separación/diversidad; conserva valores y
coincide en bytes con el blob Git. Ambos dry-runs terminan sin trabajos pendientes.

## Bloqueos y limitaciones

- DADA2 original mantiene un desempate no controlado por la semilla R. El
  workflow actual usa un adaptador explícito que reproduce soportes y tablas
  entre procesos/hilos en el entorno registrado; no se extrapola esa garantía
  al paquete original, a otros entornos ni a otro orden de entradas. Las
  asignaciones cercanas a 50 siguen siendo una sensibilidad, no evidencia primaria.
- Siete estudios permanecen pendientes por falta o inconsistencia de evidencia
  pública; no se inventan primers, tiempos, tratamientos ni asignaciones de run.
- Hay taxonomía y diversidad descriptiva del primer estudio completo; faltan
  los otros tres estudios y la inferencia biológica. Sus 93 muestras representan
  cuatro lotes, no 93 fermentaciones independientes. El piloto solapa seis de ellas.
- En el piloto, el 20,99 % de las lecturas pertenece a ASVs etiquetadas como
  orgánulos; en el estudio completo son 157.550 lecturas (6,99 %), en 18 ASVs.
  La cobertura de las tablas originales incluye esas lecturas en su denominador.
  Su separación está configurada y registrada en la tabla derivada; las
  etiquetas no identifican especie huésped ni demuestran contaminación.
- Profundidad desigual y dependencia del pseudoconteo respecto a la escala de
  conteos limitan las comparaciones descriptivas. La riqueza no está rarefactada;
  no se calcula significación ni se elige la sensibilidad por separación visual.
- Las ASVs del piloto de 237–243 nt y las completas de 236–247 nt no se filtran
  por una longitud canónica; el depósito
  tiene evidencia de recorte previo y su procesamiento exacto es incompleto.
- Los residuos raros en orientación inversa y la ausencia de controles negativos
  requieren interpretación posterior; no demuestran contaminación ni su ausencia.
- País, estudio, tecnología y región 16S pueden estar confundidos; no se les
  atribuirán efectos separables sin respaldo del diseño.
- Licencia y metadatos de autoría del repositorio siguen pendientes del titular.

## Próximas tareas

1. Registrar agregación de submuestras, tratamiento de ceros y modelos
   longitudinales antes de estimar efectos; conservar lote y sensibilidades
   temporal absoluta/relativa, sin utilizar el piloto para inferencia.
2. Procesar los otros estudios por separado, revisando primers, V3-V4 y PacBio
   antes de armonizar taxonomía e integrar efectos y heterogeneidad.

## Control de versiones

El hito piloto se desarrolló en `codex/pilot-diversity`, creada desde `main`
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
La ampliación usa `codex/full-study-processing`, creada desde `868afa8`. No se
declara el hito multistudio: solo el primer estudio tiene ejecución descriptiva
completa. El bloqueo R2 se resolvió; la continuación usa `codex/full-study-resume`.
El hito `v0.5.1-full-study` consolida esta ampliación, tablas de composición y
figuras temporales, con la limitación del RNG taxonómico documentada. La entrega
está en [PR #6](https://github.com/Arnoldd-h/Cacao-Fermentation-Microbiome/pull/6).

La continuación `codex/taxonomy-seeded-ties` registra el método en `e008f2b`
y los resultados en `a5af049`. El hito `v0.5.2-taxonomy-reproducibility` se
entrega mediante [PR #7](https://github.com/Arnoldd-h/Cacao-Fermentation-Microbiome/pull/7),
limitado a reproducción taxonómica del primer estudio y su piloto. No se
declara todavía un hito multistudio. Ambos workflows quedan sin trabajos
pendientes; la siguiente fase es registrar y contrastar el modelo longitudinal.
