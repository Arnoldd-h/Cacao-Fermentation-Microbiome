# Procesamiento completo de PRJNA492720

## Alcance registrado

`8fd5334` registra la ampliación y `b10a15e` la elegibilidad por formato,
antes de observar los resultados del estudio completo. Se auditan todos los
94 runs incluidos y se comparan accessions, layout, estrategia, alias,
URLs, bytes, MD5, número de reads y bases contra ENA, además del BioProject NCBI.
Una discrepancia obliga a revisar el inventario y detiene la generación.

El inventario conserva 94 runs. El manifest de procesamiento contiene 93
runs con streams `_1` y `_2` completos: cuatro lotes, 47 observaciones lote-tiempo
y dos estratos, sin agregar submuestras. Los tiempos son medidas repetidas;
los runs no equivalen a fermentaciones independientes. Campos desconocidos
permanecen vacíos y se conservan las fuentes originales de metadata.

La columna histórica `pilot_order` se mantiene por compatibilidad con el
lector validado; aquí representa el orden completo del manifest, no una selección
de seis muestras. Los IDs ASV son locales a cada ejecución. El aprendizaje de
errores se repite para el estudio completo; no se fusiona con la tabla piloto.

## Auditoría de empaquetado

Antes de procesar se inspeccionaron los siete archivos de tres runs con
empaquetado excepcional, verificando bytes, MD5, gzip, estructura y cabeceras:

| Run | Archivos pareados | Lecturas en stream aislado | Decisión |
|---|---:|---:|---|
| SRR7899650 | 2 | 1.486 | Usar el par y registrar el stream aislado |
| SRR7899730 | 2 | 1.252 | Usar el par y registrar el stream aislado |
| SRR7899803 | 0 | 11.064 | Excluir sólo de la ruta paired-end |

El último archivo tiene 11.064 IDs únicos: 11.056 sufijos `/1` y ocho `/2`.
No hay pares completos ni streams `_1`/`_2`; no se inventan mates ni se mezclan
ASVs inferidas single-end con ASVs ensambladas. La decisión responde al
formato público antes del análisis, sin selección por abundancia o significación.

`results/planning/pacheco_montealegre_2020_colombia/fastq_stream_log.tsv`
registra los 189 streams públicos, 186 usados y tres excluidos, con URLs,
bytes, MD5 y razones. Los 186 archivos usados suman 394.515.770 bytes;
los streams aislados suman 872.189 bytes. Las auditorías de fuentes y de
empaquetado tienen registros de procedencia independientes.

## Métodos y ejecución

`workflow/study.Snakefile` reconstruye manifest, FASTQ, QC, DADA2, taxonomía,
separación bacteriana y diversidad descriptiva, con validación en cada etapa.
`config/full_study.yaml` fija alcance, rutas y recursos operativos; no altera
primers, filtros ni decisiones de `config/config.yaml`. La reserva de disco
equivale a seis veces los bytes comprimidos y es una estimación operativa.
Se limitan cuatro tareas simultáneas y 4.500 MB; DADA2 y taxonomía usan dos hilos.

QC verifica íntegramente bytes, MD5 y FASTQ pareados antes de FastQC/MultiQC.
Cutadapt conserva la política sin descartar no coincidencias y sin trimming
de calidad. La búsqueda de primers conserva ocho casos por dirección/run,
antes y después del recorte. El cierre QC reconcilia conteos, muestras,
parámetros ejecutados y cada caso de búsqueda; una discrepancia impide DADA2.

DADA2 conserva los parámetros del [piloto](dada2_methods.md). SILVA NR99 138.2,
bootstrap 80 y sensibilidad 50 conservan la [política taxonómica](taxonomy_methods.md).
La separación y los índices alfa, CLR/Aitchison (pseudoconteos 1 y 0,5),
Bray-Curtis sobre proporciones y PCA conservan los [métodos de diversidad](diversity_methods.md).
No hay rarefacción ni pruebas temporales, PERMANOVA, abundancia diferencial o
meta-análisis. La ampliación describe las bibliotecas; aún no prueba sucesión.

Las figuras completas omiten las 93 etiquetas superpuestas. Alfa muestra
el índice de muestra en el orden tabular; no es un eje temporal. PCA mantiene
igual escala geométrica por eje y muestra ambas sensibilidades. Todas las
muestras se identifican en las tablas correspondientes; no se trazan líneas
que impliquen trayectorias entre lotes o estratos. Exportación PDF/SVG/PNG, 300 dpi.

Cada cierre científico registra commit, configuración, versiones, hashes de
entradas y salidas y estado Git. `git_dirty=true` puede reflejar los nuevos
resultados aún sin commit; ese estado se registra sin ocultarlo. Los crudos,
la base SILVA, los informes completos de herramientas y RDS quedan fuera de
Git. Los resúmenes, tablas derivadas y figuras pequeñas se versionan.

## Disponibilidad y condición de avance

`scripts/study/audit_downloads.py` verifica íntegramente los archivos locales
por bytes, MD5, gzip/CRC y estructura FASTQ. Para los archivos faltantes consulta
las cabeceras HTTP y registra la URL efectiva y el tipo de contenido. Escribe
`download_status.json` y su procedencia; salida 2 indica datos incompletos.
Una descarga fallida no añade exclusiones al manifest ni habilita DADA2 con
una selección de runs. Los datos verificados se reutilizan al reanudar.

El 2026-10-09, los 186 enlaces se inspeccionaron. La URL de
`SRR7899804_2.fastq.gz` devolvió un directorio HTML vacío, mientras la API ENA
siguió declarando 1.458.078 bytes y MD5 `49b5d1b26dd75169352b52aace2c5524`.
Se repitieron accesos HTTP/HTTPS, un segundo hostname público y FTP. La API
no ofreció archivos submitted/SRA alternativos para ese run. Esta discrepancia
de disponibilidad es distinta de la elegibilidad de SRR7899803: SRR7899804
permanece en el manifest de 93 runs y el procesamiento completo espera su R2.
Los mecanismos oficiales de transferencia se describen en la
[documentación ENA](https://ena-docs.readthedocs.io/en/latest/retrieval/file-download.html).
