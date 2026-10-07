# Piloto DADA2: método, validaciones y alcance

## Datos y objetivo

El piloto procesa exclusivamente las seis muestras de PRJNA492720 definidas en
`metadata/pilot_manifest.tsv`. Son amplicones V4 Illumina paired-end del estudio
Pacheco-Montealegre et al. (2020). Las entradas son los archivos posteriores a
Cutadapt, bajo `data/interim/{study_id}/pilot`; los archivos crudos permanecen
inmutables. Este ensayo comprueba el funcionamiento y la retención del pipeline.
No estima efectos temporales ni convierte las seis muestras en seis fermentaciones
independientes. Conserva la identidad del lote, el tiempo y la etapa de cada muestra.

## Parámetros fijados antes de observar ASVs

Los parámetros científicos se guardan en
`config/config.yaml::amplicon_processing.PRJNA492720.dada2` y su decisión en
`protocol/analysis_decisions.md`. El commit de adopción inicial es `44f5f4a`.

| Operación | Parámetros | Fundamento |
|---|---|---|
| Filtrado pareado | `truncLen=c(0,0)`, `trimLeft=c(0,0)`, `maxN=0`, `maxEE=c(2,2)`, `truncQ=2`, `minLen=150`, eliminación PhiX e identificación de pares activas | El QC del piloto muestra medianas de 222–223 nt y calidad suficiente para conservar el solapamiento. El límite de 150 nt evita conservar fragmentos muy cortos tras truncamiento por calidad; no descarta ninguna longitud observada antes del filtrado. |
| Aprendizaje de errores | Un modelo por dirección, `nbases=1e8`, sin aleatorizar archivos, `MAX_CONSIST=10`, `OMEGA_C=0` | Emplea todas las bases disponibles si el piloto no alcanza el objetivo de bases. No confunde dicho objetivo con un mínimo obligatorio para ejecutar. |
| Inferencia | Cada muestra independiente, `pool=FALSE`, `OMEGA_A=OMEGA_C=1e-40` | Mantiene muestras y estudios separados. |
| Ensamblaje | Solapamiento mínimo 12 nt, cero discrepancias, sin concatenación y sin recorte de overhang | Exige respaldo de ambas lecturas para cada secuencia ensamblada. |
| Quimeras | Consenso dentro del estudio, `minFoldParentOverAbundance=1.5`, `minParentAbundance=2`, `minSampleFraction=0.9`, `ignoreNNegatives=1`, sin variantes a una diferencia | Parámetros por defecto verificados en DADA2 1.34.0 para `isBimeraDenovoTable`. |
| Longitud ensamblada | Sin filtro adicional | Se exporta la distribución observada. Un filtro posterior requiere evidencia y una nueva decisión explícita. |

FastQC muestra medias mínimas por intervalo de posiciones de 28,06–34,36 para R1
y 31,63–32,65 para R2; las medianas mínimas son ≥31. No se observó un descenso
abrupto de calidad que justificara truncar todas las lecturas a una posición fija.
Las siete lecturas R2 recortadas por Cutadapt amplían el rango observado a 185–223
nt; las medianas siguen siendo 222–223 nt. Para un amplicón V4 de aproximadamente
250 nt, dos fragmentos de 150 nt permiten aproximadamente 50 nt de solapamiento;
esta expectativa orienta el guardarraíl y no sustituye al ensamblaje real.

El [tutorial oficial de DADA2](https://benjjneb.github.io/dada2/tutorial.html)
respalda el filtrado por errores esperados, el aprendizaje separado de errores y
el seguimiento de lecturas a través de inferencia, ensamblaje y quimeras. El
[manual oficial](https://www.bioconductor.org/packages/release/bioc/manuals/dada2/man/dada2.pdf)
documenta los argumentos de estas funciones; las firmas de la instalación 1.34.0
se verificaron localmente antes de fijarlos. La versión del manual en línea puede
cambiar; la ejecución conserva versiones y parámetros efectivos.

## Ejecución reproducible

Desde Linux/WSL, en la raíz del repositorio:

```bash
python3 scripts/environment/run_in_environment.py snakemake --snakefile workflow/Snakefile --cores 2 pilot_dada2
```

El target reconstruye QC y ejecuta un validador independiente que comprueba
conteos, secuencias, longitudes y checksums de salidas y entradas originales.
Para invocar únicamente la etapa R con entradas QC existentes:

```bash
python3 scripts/environment/run_in_environment.py Rscript scripts/dada2/run_pilot.R \
  --config config/config.yaml \
  --manifest metadata/pilot_manifest.tsv \
  --output-dir results/dada2/pilot \
  --intermediate-dir results/intermediate/dada2/pilot \
  --raw-quality results/qc/pilot/raw_read_quality.tsv \
  --trimmed-quality results/qc/pilot/trimmed_read_quality.tsv \
  --threads 2
```

Agregar `--validate-only` valida metadata, parámetros y existencia de las entradas
sin filtrarlas ni generar resultados. Las pruebas sintéticas de contratos se
ejecutan mediante:

```bash
python3 scripts/environment/run_in_environment.py Rscript tests/dada2/test_helpers.R
python3 scripts/environment/run_in_environment.py Rscript tests/dada2/test_filter.R
```

El manifiesto debe contener un único estudio/BioProject, lecturas paired-end
Illumina e identificadores únicos. El programa rechaza conteos inconsistentes,
metadata temporal inválida y muestras sin lecturas filtradas, ensambladas o sin
quimeras. No excluye estas muestras silenciosamente. Los warnings permanecen
visibles y se archivan en `warnings.txt`; una ejecución fallida no recibe el
marcador `SUCCESS`.

El argumento `--threads` puede reducir los hilos hasta uno, sin superar el límite
configurado; la procedencia distingue los parámetros del estudio del número
efectivo de hilos. `error_learning.tsv` conserva las rondas del aprendizaje, la
repetición de matrices de error y si se alcanzó `MAX_CONSIST`. Se reproduce el
criterio de parada de DADA2 1.34: una matriz repetida antes del límite, incluido
un posible ciclo. Tasas de error positivas por sí solas no demuestran convergencia.
Una corrida que termine con warnings o sin convergencia requiere revisión aunque
haya generado tablas.

## Productos y trazabilidad

- `read_tracking.tsv`: entrada cruda si se aportan tablas QC, entrada tras Cutadapt,
  filtrado, denoising por dirección, ensamblaje, quimeras, proporciones y ASVs por
  muestra. Comprueba conservación de conteos en cada transición.
- `asv_counts.tsv`, `asv_sequences.tsv` y `asv_sequences.fasta`: conteos enteros y
  secuencias sin asignaciones taxonómicas inventadas. Los IDs contienen el estudio
  y se ordenan por secuencia; son identificadores de esta ejecución, no claves
  globales para integrar estudios o distintas versiones de la tabla.
- `sequence_length_distribution.tsv` y `merge_diagnostics.tsv`: longitudes antes
  y después de quimeras y pares aceptados/rechazados durante el ensamblaje.
- `error_model_F/R.pdf`, `.svg` y `.png`: curvas diagnósticas de errores observados
  y estimados; PNG a 300 dpi. Son diagnósticos técnicos, no figuras biológicas de
  publicación.
- `run_provenance.yaml`, `config_snapshot.yaml`, `input_checksums.tsv`,
  `output_checksums.tsv`, `software_versions.tsv` y `session_info.txt`: fecha UTC, commit, estado Git,
  comando, parámetros efectivos, semilla y versiones. Los MD5 de configuración,
  manifiesto, scripts, QC y FASTQ se calculan antes y se comprueban al terminar.
- `results/intermediate/dada2/pilot/`: FASTQ filtrados y objetos RDS regenerables,
  excluidos de Git.

Los resultados del piloto no sustituyen controles negativos ausentes, no asignan
taxonomía, no prueban ausencia de contaminación y no validan de forma automática
los parámetros para otro estudio, plataforma o región 16S. Antes de escalar se
revisan retención, ensamblaje, longitudes y ajuste de errores; cualquier cambio
metodológico exige configuración, decisión registrada y commit específico.

## Corrección del identificador SRA

El primer intento del 4 de octubre de 2026 falló antes de filtrar por no poder
autodetectar los identificadores de pares SRA. Las cabeceras reales tienen la
forma `@SRR...1 1/1` y `@SRR...1 1/2`; el primer campo identifica al mismo par.
Se configuraron `id.field=1` e `id.sep="\\s"`, manteniendo `matchIDs=TRUE`. La
corrección quedó en el commit `eafadd7` y tiene una regresión con FASTQ sintéticos
pareados. No se modificaron umbrales científicos como respuesta al fallo.
El intento fallido se conserva en
`results/dada2/attempts/sra_identifier_failure.yaml`.

## Resultado técnico observado

La primera ejecución técnica completa del 7 de octubre de 2026 (UTC), previa a
la consolidación del código y la ejecución final del workflow, conservó 165.034 de
176.798 pares durante el filtrado (93,35 %), ensambló 160.649 y retuvo 156.561
pares sin quimeras (88,55 % de la entrada). Las 131 secuencias ensambladas se
redujeron a 93 ASVs después del consenso de quimeras. Ninguna de las seis muestras
quedó vacía; la retención final varió de 80,01 % a 94,39 %.

Ambos modelos de error convergieron en cinco rondas, por debajo del máximo de
diez. La revisión visual de ambas curvas mostró que los ajustes siguen la
tendencia de las tasas observadas, con ruido en transiciones infrecuentes y
diferencias respecto a las tasas nominales de Phred. Se conservaron doce warnings
de representación `log-10 transformation introduced infinite values`, debidos a
frecuencias observadas nulas; no son avisos de fallo de convergencia. La existencia
de curvas y convergencia no equivale a una validación biológica del modelo.

Las ASVs retenidas tienen 237–243 nt; esta distribución estrecha se conserva sin
forzarla a una longitud V4 canónica. Los FASTQ depositados ya tenían evidencia de
procesamiento previo, cuya secuencia exacta no puede reconstruirse por completo
desde este piloto. La pertenencia taxonómica y la posible contaminación siguen
pendientes de clasificación y controles apropiados. Los mensajes de secuencias
duplicadas se emitieron durante el ensamblaje con `returnRejects=TRUE`, que conserva
también los emparejamientos rechazados para el diagnóstico. Sólo los registros con
`accept=TRUE` se pasan a `makeSequenceTable`; la validación independiente comprueba
que las abundancias agregadas conservan los conteos.

Las cifras proceden de `results/dada2/pilot/summary.tsv`, `read_tracking.tsv`,
`sequence_length_distribution.tsv` y `error_learning.tsv`. La ejecución definitiva
del workflow conserva su propio commit y hashes en `run_provenance.yaml`.
