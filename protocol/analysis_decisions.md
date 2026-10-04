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

## Pendiente

### Identificadores pareados SRA, 2026-10-03

El primer intento de `filterAndTrim` detectó que la identificación automática
no reconocía encabezados como `@SRR7899687.1 1/1` y `@SRR7899687.1 1/2`.
Se configura `id.field=1` e `id.sep="\\s"`: el primer token contiene el
identificador compartido de la pareja. `matchIDs=true` se mantiene y comprueba
la correspondencia real; no se permite desactivarlo para sortear el error.
Esto corrige parsing del depósito, sin modificar filtros, umbrales ni lecturas.
El caso se incorpora a las pruebas con encabezados sintéticos del mismo formato.

### Decisiones aún abiertas

- Versión de SILVA y método definitivo de clasificación.
- Validar los parámetros iniciales de DADA2 con los diagnósticos del piloto antes de escalar.
- Tratamiento de ceros y pseudoconteo para CLR.
- Método definitivo de abundancia diferencial según el diseño disponible.
- Selección de licencia del repositorio.
