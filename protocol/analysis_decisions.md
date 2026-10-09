# Decisiones de análisis

| Fecha | Decisión | Justificación | Alternativas consideradas | Impacto esperado | Estado |
|---|---|---|---|---|---|
| 2026-08-19 | Limitar la primera fase a amplicones bacterianos 16S rRNA. | Evita mezclar tecnologías y permite estabilizar un pipeline defendible. | Integrar ITS o shotgun desde el inicio. | Menor alcance inicial, mayor comparabilidad. | Adoptada |
| 2026-08-19 | Procesar ASVs por estudio y armonizar después a género o familia. | Primers, regiones y protocolos diferentes no producen ASVs directamente comparables. | Concatenar tablas ASV cross-study. | Reduce falsa equivalencia y conserva efectos dentro de estudio. | Adoptada |
| 2026-08-19 | Mantener tiempo absoluto y relativo; derivar etapas desde configuración. | Las fermentaciones tienen duraciones distintas y los límites deben someterse a sensibilidad. | Hardcodear `early/mid/late` en scripts. | Facilita modelos continuos y análisis de sensibilidad. | Adoptada |
| 2026-08-19 | Usar PRJNA492720 como candidato piloto tras verificar su subconjunto publicado de 94 runs 16S. | Tiene diseño longitudinal espontáneo, V4, paired-end y metadata temporal recuperable. | Empezar con todos los proyectos o con WGS. | Permite una vertical slice contenida antes de descargas masivas. | Adoptada |
| 2026-08-19 | Mantener PRJNA865318 y PRJNA1104253 como pendientes. | El primero tiene tiempo incompleto/ambiguo; el segundo mezcla diseños y su subconjunto 16S verificable es controlado in vitro. | Incluirlos automáticamente por contener `AMPLICON`. | Evita ampliar el análisis primario sin evidencia suficiente. | Adoptada |
| 2026-08-19 | Excluir PRJNA552479, PRJNA1257864 y PRJNA1264670 de Fase I. | ENA reporta WGS, no amplicones bacterianos 16S. | Procesarlos como si fueran 16S o ampliar ahora a shotgun. | Mantiene el alcance y deja trazabilidad para fases futuras. | Adoptada |
| 2026-08-19 | Usar una consulta NCBI BioProject documentada y ENA para descubrir candidatos. | Una lista manual de accessions no demuestra cobertura sistemática ni es fácilmente actualizable. | Mantener solo seis estudios semilla. | Produce un universo reproducible de 34 proyectos y una cola auditable de revisión. | Adoptada |
| 2026-08-19 | Impedir que la triaje automática promueva estudios a inclusión. | Las señales textuales no verifican primers, temporalidad, espontaneidad ni coherencia con la publicación. | Incluir automáticamente todo proyecto `AMPLICON`. | Reduce falsos positivos; conserva 18 candidatos como revisión manual. | Adoptada |
| 2026-08-19 | Incluir PRJNA627078 como segundo estudio primario. | El artículo y ENA verifican fermentación espontánea en cuatro cajas, 16S V3-V4 separable, primers publicados, cinco tiempos de 0 a 120 h y FASTQ paired-end. | Mantenerlo pendiente por coexistir con 26S y controles de cosecha limpia. | Añade 60 corridas tradicionales; 26S y seis controles limpios quedan excluidos explícitamente. | Adoptada |
| 2026-08-19 | Mantener PRJEB82327 y PRJEB82871 pendientes pese a identificar sus controles espontáneos full-length 16S. | ENA permite aislar 48 y 16 corridas PacBio longitudinales, respectivamente, pero aún no hay artículo vinculado ni primers verificados. | Incluirlos usando solo la descripción del depósito. | Conserva los subconjuntos candidatos sin incorporarlos como evidencia primaria. | Adoptada |
| 2026-08-19 | Mantener PRJEB53853 pendiente y excluir PRJNA962540. | En PRJEB53853 falta mapear marcador y tratamiento por run; PRJNA962540 concentra siete amplicones en un único BioSample sin tiempo y solo tiene dos tiempos WGS. | Inferir marker/tratamiento o tratar WGS como 16S. | Evita pseudorreplicación e inclusión de tecnologías incompatibles. | Adoptada |
| 2026-08-19 | Seleccionar PRJNA492720 como piloto y preparar seis muestras sin iniciar aún la descarga. | La comparación configurada lo sitúa primero: metadata alta, primers y paper verificados, paired-end, 13 tiempos, cuatro fermentaciones y 377,07 MiB para el subconjunto completo. | PRJNA627078, también elegible, tiene cinco tiempos y 1.320,17 MiB; los pendientes no reciben ranking. | El manifest cubre dos lotes independientes por etapa y exige dos FASTQ, bytes y MD5 declarados por run. | Adoptada |
| 2026-08-19 | Descargar únicamente los seis runs del manifest piloto y exigir bytes, MD5, gzip y estructura FASTQ válidos. | La vertical slice suma 25.154.687 bytes, existe espacio suficiente y `data/raw/` está excluido de Git. | Esperar al inventario completo o descargar las 94 corridas. | Habilita QC real con riesgo y almacenamiento contenidos; 12/12 archivos fueron validados. | Adoptada |
| 2026-08-20 | Ejecutar el stack Bioconda en Ubuntu/WSL2 con prioridad estricta y fijar dependencias con ABI sensible. | Bioconda no soporta Windows nativo; las resoluciones sin restricciones permitían instalar combinaciones que fallaban al cargar DADA2 o ANCOMBC. | Instalar herramientas heterogéneas en Windows, aceptar solo que el solver termine o usar versiones no fijadas. | Hace reproducibles FastQC, Cutadapt, Snakemake, R/DADA2 y los métodos posteriores; las restricciones se reevaluarán al actualizar DADA2 o ANCOMBC. | Adoptada |
| 2026-08-25 | Conservar completas las lecturas de 214-223 nt durante la verificación de primers, sin fijar todavía truncamiento DADA2. | Los 12 FASTQ crudos pasan calidad por base, tienen GC de 53-56 %, cero reads marcados como mala calidad y el módulo de adaptadores pasa; las fallas de composición, duplicación y sobrerrepresentación son compatibles con bibliotecas de amplicones. | Truncar por costumbre antes de comprobar primers y solapamiento. | Mantiene la máxima longitud disponible para verificar orientación y calcular el solapamiento real antes de DADA2. | Adoptada |
| 2026-08-25 | Configurar para PRJNA492720 los núcleos 515F `GTGCCAGCMGCCGCGGTAA` y 806R `GGACTACHVGGGTWTCTAAT`, junto con sus constructos L1/L2 completos de Table S2. | El artículo primario identifica V4 y primers modificados 515F-806R; `Data_Sheet_2.PDF`, Table S2, resuelve que la modificación son linkers 5′ y proporciona las secuencias exactas. | Usar las variantes posteriores Parada/Apprill solo por compartir los nombres 515F/806R o conservar descripciones no secuenciales. | Evita recortar una variante de primer incorrecta y deja fuente primaria y constructos auditables en configuración. | Adoptada |
| 2026-08-25 | Tratar los FASTQ depositados del piloto como lecturas ya recortadas de primers y conservar las lecturas sin coincidencia en Cutadapt. | En 176.798 reads por dirección no se detectó 515F en R1; solo 7 coincidencias de 806R en R2 (0,004 %) y ningún constructo L1/L2 completo en la orientación esperada. Se observaron residuos raros en orientación inversa, principalmente en SRR7899884, coherentes con artefactos/adaptadores de una fracción menor y no con primers conservados de forma sistemática. El artículo además declara trimming con Trimmomatic en su procesamiento publicado. | Descartar todas las lecturas sin primer o aplicar trimming obligatorio como si los FASTQ fueran crudos de instrumento. | El paso Cutadapt será no destructivo para lecturas sin coincidencia y deberá demostrar retención, sin asumir que el primer sigue presente. | Adoptada |
| 2026-08-25 | Ejecutar Cutadapt 5.0 solo para los núcleos verificados en la orientación esperada, con coincidencia de longitud completa, `error_rate=0.1`, sin indels, sin quality trimming y sin descartar untrimmed. | La detección previa mostró primers sistemáticamente ausentes; mezclar recorte de calidad o descarte impediría atribuir las pérdidas. La ejecución conservó 176.798/176.798 pares, recortó solo 7 R2 y retiró 238 bases. | Descartar untrimmed, recortar orientaciones raras o combinar primers y calidad en un paso. | Produce FASTQ interim pareados y una tabla por dirección; cualquier pérdida posterior podrá atribuirse a QC/DADA2 y no a descarte de Cutadapt. | Adoptada |
| 2026-10-03 | Fijar parámetros iniciales DADA2 1.34.0 del piloto antes de ejecutarlo: sin truncamiento fijo, maxN=0, maxEE=(2,2), truncQ=2, minLen=150, aprendizaje por dirección hasta 10 iteraciones y 100 millones de bases disponibles, inferencia independiente y merging con 12 nt de solapamiento sin discrepancias. | Los FASTQ post-Cutadapt miden 185–223 nt; los perfiles por ciclo no muestran un colapso abrupto de calidad y ya se verificó recorte previo al depósito. minLen=150 es un límite técnico inicial para lecturas excepcionalmente acortadas por truncQ; no descarta las longitudes actualmente observadas. Los valores restantes siguen el tutorial/manual oficial, cotejados con la versión instalada. | Truncar a una longitud arbitraria, relajar merging para elevar retención o ajustar parámetros según taxones/efectos. | Se registran retención por etapa, errores aprendidos, rechazos de merging y longitudes; no se interpreta el piloto como evidencia de sucesión ni se escala antes de revisar estos diagnósticos. Semilla 20260819; no se mezclan estudios. | Adoptada para piloto técnico |
| 2026-10-03 | Usar eliminación de quimeras consensus con minFoldParentOverAbundance=1.5, minParentAbundance=2, minSampleFraction=0.9, ignoreNNegatives=1 y allowOneOff=false; no aplicar todavía filtro por longitud del amplicón ensamblado. | Son parámetros explícitos de la versión instalada; las longitudes ensambladas reales aún no se conocen. OMEGA_C=0 corresponde al aprendizaje de errores y OMEGA_C=1e-40 a la inferencia final. | Confundir los parámetros de aprendizaje e inferencia o imponer un intervalo post-merge sin observar los datos. | Cualquier filtro posterior exigirá evidencia, configuración, validación y decisión específica, conservando las tablas sin filtrar. | Adoptada para piloto técnico |

### Estructura del diseño, 2026-10-03

Se fija en `config.analysis_design` que la unidad independiente es el lote;
los tiempos son medidas repetidas y las extracciones/estratos son submuestras
anidadas. La alternativa de contar cada run como réplica independiente se
rechaza por pseudorreplicación. Se exportará el mapa estudio/lote/hora sin
alterar abundancias y se conservará la identificación de cada submuestra para
la agregación y sensibilidad posteriores. Esta decisión se basa en el diseño
documentado en `reports/dataset_inventory.md`, no en resultados de abundancia.
El piloto técnico no habilita inferencia biológica. Los moderadores asociados
con estudio/región 16S se declararán no identificables mientras el diseño no
permita separarlos; horas absolutas y tiempo relativo se conservan para evaluar
la dependencia respecto del último tiempo observado de cada lote.

Fuentes de DADA2: [tutorial oficial](https://benjjneb.github.io/dada2/tutorial.html)
y [manual de referencia](https://www.bioconductor.org/packages/release/bioc/manuals/dada2/man/dada2.pdf).
Los valores por defecto se contrastaron con las funciones de DADA2 1.34.0
instaladas; la versión de referencia en línea puede ser posterior.

### Identificadores pareados SRA, 2026-10-03

El primer intento de `filterAndTrim` detectó que la identificación automática
no reconocía encabezados como `@SRR7899687.1 1/1` y `@SRR7899687.1 1/2`.
Se configura `id.field=1` e `id.sep="\\s"`: el primer token contiene el
identificador compartido de la pareja. `matchIDs=true` se mantiene y comprueba
la correspondencia real; no se permite desactivarlo para sortear el error.
Esto corrige parsing del depósito, sin modificar filtros, umbrales ni lecturas.
El caso se incorpora a las pruebas con encabezados sintéticos del mismo formato.

### Revisión completa del inventario, 2026-10-06

La revisión pública iniciada el 4 de octubre y terminada el 7 de octubre UTC
(6 de octubre en Colombia) configura los 34 candidatos: cuatro incluidos,
siete pendientes y 23 excluidos. Las fuentes y alternativas por candidato
quedan en `config/datasets.yaml` y `reports/inventory_review_2026-10-04.md`.

- Se incluyen los controles espontáneos F1/F2 de PRJEB40850 (12 runs, V4,
  0–92 h) y F01/F02 de PRJEB57747 (16 runs, 16S completo, 0–120 h). Los
  métodos primarios y archivos depositados con marcador explícito permiten
  distinguirlos de inoculados, ITS y WGS. Se rechazó seleccionar por tamaño
  de archivo o por la etiqueta genérica AMPLICON. Los primers proceden de los
  protocolos citados y se preservan con sus variantes exactas en configuración.
- PRJEB57747 requiere procesamiento PacBio específico; su inclusión científica
  no autoriza reutilizar parámetros Illumina. Los experimentos Costa Rica
  2017/2019 se mantienen separados y se evaluará dependencia por sitio/laboratorio.
- Las decisiones pendientes iniciales de PRJNA865318 y PRJNA1104253 quedan
  sustituidas por exclusión: sólo hay día 1 bacteriano verificable en el primero;
  los amplicones del segundo pertenecen a ensayos controlados y el componente
  natural es WGS. Se rechazó inferir FASTQ de filas de muestreo planificadas.
- PRJEB53853 permanece pendiente: se resolvieron tratamientos F1/F2, pero no la
  separación bacteriana/fúngica de bibliotecas mixtas. PRJNA420946 conserva los
  registros a 144 h, incompatibles con las 120 h de la publicación; no se
  recodifican ni se excluyen para obtener concordancia. Los otros cinco
  pendientes conservan el requisito preciso de evidencia en el informe.
- Un lote desconocido queda vacío; `study_id` ya no lo sustituye. Se preservan
  horas absolutas conocidas sin fabricar duración, tiempo relativo o etapa.
  Los valores `unknown`/`not reported` se trasladan a notas y el campo observado
  queda vacío. Esta corrección evita crear réplicas o covariables sin evidencia.

Las reglas se aplican a 2.357 runs y conservan 182 incluidos, 93 pendientes y
2.082 excluidos, con 12 lotes entre los incluidos. La selección prerregistrada
mantiene PRJNA492720 y el mismo manifest de seis runs. Se validan las tablas,
selectores de marcador/plataforma, ausencia de lotes inventados y conservación
de los tiempos en conflicto antes de consolidar el cambio metodológico.

## Decisiones aún abiertas

- Métodos posteriores de agregación taxonómica entre estudios.
- Revisar el comportamiento de DADA2 al escalar el estudio completo.
- Tratamiento de ceros y pseudoconteo para CLR.
- Evaluar sustitución de ceros para inferencia al ampliar el estudio; la política
  técnica del piloto se fija abajo y no establece el método inferencial definitivo.
- Método definitivo de abundancia diferencial según el diseño disponible.
- Selección de licencia del repositorio.

### Clasificación taxonómica del piloto, 2026-10-09

Se fija antes de observar las asignaciones SILVA NR99 138.2, archivo oficial
toGenus para DADA2 del [depósito 14169026](https://doi.org/10.5281/zenodo.14169026).
Su MD5 de origen es `1764e2a36b4500ccb1c7d5261948a414` y tamaño 139.996.892 bytes;
se comprobarán gzip, formato FASTA y SHA-256 local. La base pesada permanece
ignorada. Sus archivos fueron formateados con DADA2 1.35.4; se utilizará el
clasificador instalado 1.34.0 y se verificará la compatibilidad mediante ejecución.
Se usa el enlace directo del mismo archivo en Zenodo tras observar una descarga
anormalmente lenta por el endpoint API. Se conserva DOI, tamaño y MD5; el enlace
directo admite HTTP Range para reanudar la descarga del mismo archivo fijado.

`config/taxonomy.yaml` registra `assignTaxonomy`, seis rangos originales Kingdom
a Genus, comprobación de complemento inverso y semilla 20260819. Se conserva
la llamada sin máscara `raw_min_boot=0` y el soporte de cada rango. Se usa el
generador Mersenne-Twister /
Inversion / Rejection, ahora enlazados explícitamente desde la configuración.
La máscara primaria exige bootstrap 80 y se compara descriptivamente con 50
usando el mismo ajuste.
Se elige 80 para conservar soporte conservador al preparar integración a género;
50 evalúa pérdida de resolución en las secuencias cortas de 237–243 nt. Son
umbrales documentados por la [guía oficial](https://benjjneb.github.io/dada2/assign.html),
sin selección según taxones esperados ni efectos temporales. Un rango inferior
no se acepta si un ancestro carece de soporte. Se preservan los nombres SILVA,
incluidos grupos sin nombre de género válido; no se inventan sinónimos o especies.

No se asignan especies en esta fase por la resolución limitada de V4. Se
marcan etiquetas explícitas Chloroplast/Mitochondria, dominios distintos de
Bacteria y dominios desconocidos para revisión, conservando todas las ASVs y
lecturas. Esta referencia optimizada para Bacteria/Archaea no permite descartar
eucariotas por una ausencia de asignación. Cualquier exclusión posterior exigirá
una decisión y registro específicos. Los resúmenes son técnicos del piloto;
no habilitan inferencia biológica ni concatenación de ASVs entre estudios.

### Separación bacteriana y diversidad descriptiva, 2026-10-09

La clasificación primaria ya permitió observar etiquetas de orgánulos. Antes
de calcular diversidad se fija `config/diversity.yaml`: excluir de una tabla
derivada las ASVs con etiquetas exactas Chloroplast/Mitochondria en cualquiera
de los seis rangos de la llamada primaria (bootstrap 80), y Kingdom conocido
distinto de Bacteria. La razón es el alcance bacteriano del estudio; no se usa
abundancia, prevalencia, etapa ni significación. Las asignaciones de Kingdom
desconocidas se conservan y marcan; una falta de resolución no justifica
inventar una identidad. No se excluyen muestras: una muestra vacía causa fallo.

La coincidencia por rango conserva la nomenclatura SILVA y evita seleccionar
por fragmentos de nombres. Se conserva todo el resultado DADA2/taxonómico
original, con un registro por ASV y balances de lecturas por muestra. El
[ejemplo oficial de QIIME 2](https://docs.qiime2.org/2024.10/tutorials/filtering/)
documenta separación taxonómica de orgánulos; aquí no se instala QIIME ni se
aplica su búsqueda por subcadenas. La separación no es un análisis de
contaminación y no identifica especie huésped.

La diversidad se calcula sólo dentro del estudio a nivel ASV. No se agregan
submuestras ni se interpreta cada run como un lote independiente. Riqueza
observada, Shannon con logaritmos naturales, Gini-Simpson (1−Σp²), inverso de
Simpson y exp(Shannon) describen las bibliotecas retenidas sin rarefacción;
su comparación inferencial requiere tratar la profundidad de muestreo.
Las definiciones se contrastan con el [manual de vegan](https://vegandevs.github.io/vegan/reference/diversity.html).

Para el diagnóstico composicional se suma 1 a todos los conteos, se calcula
CLR y su distancia euclídea; la sensibilidad usa 0,5 sobre las mismas ASVs.
Es una sustitución técnica explícita, dependiente de escala/profundidad, no
una estimación de ceros estructurales. Ambos valores se fijan antes de observar
distancias, sin escoger el que favorezca separación temporal. La
[definición oficial de CLR](https://scikit.bio/docs/dev/generated/skbio.stats.composition.clr.html)
usa logaritmos de componentes positivos respecto de su media geométrica.
Bray-Curtis sobre proporciones se conserva como diagnóstico secundario.

PCA de CLR usa centrado por ASV y no estandariza varianzas; se guardan todas
las componentes hasta n−1, su varianza y las distancias originales. Las
figuras muestran muestras individuales, sin unir lotes/estratos en una
trayectoria inventada. Semilla 20260819 registrada; las operaciones actuales
son deterministas. No se ejecutan PERMANOVA, pruebas de etapas, abundancia
diferencial, core ni meta-análisis en esta rebanada técnica de seis runs.
