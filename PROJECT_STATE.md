# Estado del proyecto

Última actualización: 2026-10-31.

## Fase actual

Segundo estudio (PRJNA627078, México) adquirido y con QC basal de primers 5'
completo: 60 bibliotecas, cuatro cajas, cinco tiempos y 20 combinaciones
caja/tiempo. Los 120 FASTQ suman 1.379.724.114 bytes y pasan bytes/MD5,
gzip/CRC/estructura. Se conservan 5.689.175 pares y todas las bibliotecas;
552 pasos completados, 340 entradas/11 artefactos validados desde `9f0fd83`.

El QC detecta una limitación que debe resolverse antes de DADA2: calidad por
base pasa de 101 pass/11 warn/8 fail a 58/1/61, incluidos los 60 R2 con fail;
los 120 archivos mantienen pass por adaptadores. Persisten coincidencias con
el primer opuesto: 101.479 R1 y 83.624 R2. No se excluyen muestras ni se asigna
identidad taxonómica a esos fragmentos.

El diagnóstico `75407e9`, sin modificar lecturas, selecciona por mate el mayor
recuento residual. En SRR11579065, 12.400/13.892 coincidencias R1 y
10.342/12.295 R2 son terminales (89,26/84,12 %); no son estimaciones de toda
la cohorte. Las longitudes variables cambian los subconjuntos que contribuyen
a las posiciones finales de FastQC y deben revisarse junto con la calidad.
DADA2 permanece candidato bajo `processing_ceiling=qc`: primero se registrará
y validará el tratamiento 3' y de fragmentos cortos. Informe:
`reports/prjna627078_qc_2026-10-31.md`.

La menor profundidad es SRR11579032/16S_cocoa_B59M: 877 pares, Criollo, caja 2,
120 h, réplica 2. Se conserva. Las etiquetas de réplica no se convierten en
estratos espaciales: `sampling_stratum` sigue vacío y la caja es la unidad
independiente. NCBI/ENA verifican también título/descripción y equipo, además
de identificadores, diseño de lectura y archivos. El inventario no cambia.

El primer estudio (PRJNA492720) conserva su ejecución completa de 93 muestras:
782 ASVs/2.254.530 pares sin quimeras; tras separar orgánulos, 764 ASVs y
2.096.980 lecturas. Su modelo longitudinal exploratorio tiene cuatro lotes,
46 tiempos/92 muestras principales y sensibilidad 47/93: cosenos de pendientes
0,249–0,584, R² 0,118–0,164 y p=0,125/q BY=0,30625 en seis escenarios.
La resolución del contraste impide rechazar a 0,05; no demuestra ausencia de
sucesión ni conservación entre estudios. Método/resultados en los informes
longitudinales y `v0.5.3-longitudinal`.

El adaptador taxonómico con semilla conserva su validación por repetición
entre procesos/hilos. Los resultados históricos mantienen sus versiones y
hashes productores; la ampliación de código compartido requiere renovar las
etapas afectadas si se reprocesan con este checkout. Sus tablas científicas
no se han editado ni se combinan ASVs entre regiones 16S.

Desarrollo: suite de 188 pruebas aprobada (179,502 s), más las regresiones
nuevas de diseño ENA, conexiones detenidas y coordenadas; 193 pruebas distintas
verificadas mediante suite y pruebas específicas. FASTQ, intermedios y reportes
pesados quedan fuera de Git. Los otros dos estudios todavía no tienen
procesamiento de lecturas; la integración/core/meta-análisis siguen pendientes.

## Completado

- Segundo estudio: manifest de 60 runs, 120 FASTQ íntegros, QC basal completo y
  diagnóstico reproducible de residuos terminales/internos, sin exclusiones.
- Descargador con timeout HTTP configurado y prueba de transferencia local
  detenida/reanudada mediante Range, conservando los archivos crudos válidos.


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
- Hay taxonomía, diversidad descriptiva y un modelo longitudinal exploratorio
  del primer estudio; el segundo sólo tiene QC basal y requiere corregir el
  preprocesamiento antes de ASVs. Faltan los otros dos y la inferencia cross-study. Sus
  93 muestras representan cuatro lotes, no 93 fermentaciones independientes.
  El contraste fijado no puede bajar de p=0,125 y supone simetría central;
  el piloto solapa seis muestras y no se usa como evidencia independiente.
- En el piloto, el 20,99 % de las lecturas pertenece a ASVs etiquetadas como
  orgánulos; en el estudio completo son 157.550 lecturas (6,99 %), en 18 ASVs.
  La cobertura de las tablas originales incluye esas lecturas en su denominador.
  Su separación está configurada y registrada en la tabla derivada; las
  etiquetas no identifican especie huésped ni demuestran contaminación.
- Profundidad desigual y dependencia del pseudoconteo respecto a la escala de
  conteos limitan las comparaciones descriptivas. La riqueza no está rarefactada;
  no incluye pruebas alfa por muestra ni selección de sensibilidad por separación
  visual. Los contrastes globales CLR tienen su propio plan por lote y límites.
- Las ASVs del piloto de 237–243 nt y las completas de 236–247 nt no se filtran
  por una longitud canónica; el depósito
  tiene evidencia de recorte previo y su procesamiento exacto es incompleto.
- Los residuos raros en orientación inversa y la ausencia de controles negativos
  requieren interpretación posterior; no demuestran contaminación ni su ausencia.
- País, estudio, tecnología y región 16S pueden estar confundidos; no se les
  atribuirán efectos separables sin respaldo del diseño.
- Licencia y metadatos de autoría del repositorio siguen pendientes del titular.

## Próximas tareas

1. Registrar y validar el tratamiento de primers 3'/fragmentos cortos y revisar
   calidad por posición en PRJNA627078 antes de DADA2, taxonomía y diversidad.
   Procesar después los dos estudios restantes con métodos apropiados a su
   región/plataforma, incluido PacBio.
2. Registrar armonización género/familia, modelos por taxón/etapa y core antes
   de la inferencia cross-study. Las pendientes ASV actuales no son combinables
   directamente con ASVs de otros estudios.
3. Evaluar si un nuevo plan temporal debe representar curvatura/transitorios,
   sin cambiar el contraste actual para buscar significación.

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
declara todavía un hito multistudio. En ese hito, ambos workflows quedaron sin trabajos
pendientes; la siguiente fase fue registrar y contrastar el modelo longitudinal.

El hito `v0.5.3-longitudinal` entrega ese primer modelo exploratorio mediante
[PR #8](https://github.com/Arnoldd-h/Cacao-Fermentation-Microbiome/pull/8):
plan `c97345c`, implementación `7b6d5a6`, renovación operativa `6415bfd` y
resultados `0132f7d`. Al cerrar ese hito, sus tres workflows quedaron sin pendientes; comenzó después
el procesamiento independiente de los otros estudios.


La continuación `codex/prjna627078-qc` registra el perfil en `d3f6149`,
la verificación de diseño/equipo en `3faac65`/`9f319c1`, el timeout en
`49c4161` y la adquisición íntegra en `9f0fd83`. El QC basal procede de
`9f0fd83` y el diagnóstico de `75407e9`. El hito
`v0.5.4-second-study-qc` conserva este control basal con sus limitaciones;
no declara una segunda ejecución DADA2 ni el hito multistudio. La renovación
de fuentes conserva el manifest (31 entradas/seis artefactos validados),
y el workflow mexicano queda sin trabajos pendientes. Los hashes de los
inputs pequeños/outputs QC y de los dos diagnósticos se comprueban al cerrar.

Entrega del cierre de datos `376bb9e` mediante
[PR #9](https://github.com/Arnoldd-h/Cacao-Fermentation-Microbiome/pull/9).
El hito conserva el QC basal y el diagnóstico de posiciones; la siguiente
unidad de desarrollo debe registrar y validar la limpieza 3' y revisar la
calidad por posición antes de habilitar DADA2.
