# Reanudación del estudio completo — 2026-10-31

## Recuperación de la fuente pública

El R2 de `SRR7899804` vuelve a estar disponible en la URL ENA registrada:

<https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR789/004/SRR7899804/SRR7899804_2.fastq.gz>.

El acceso HTTPS devuelve HTTP 200, `application/x-gzip` y 1.458.078 bytes.
La API ENA conserva las mismas URLs, tamaños y MD5 de ambos mates. El R2
descargado coincide con el MD5 original `49b5d1b26dd75169352b52aace2c5524`;
no se sustituyó la fuente ni se reconstruyeron lecturas.

El `.part` anterior contenía 797 bytes de HTML, SHA-256
`161e2e8a066b0e07fdc925c15d25adb3b169bc16186799525470c98c3520eff6`.
Se trasladó a `results/tmp/ena_html_2026-10-09_SRR7899804_2.fastq.gz.part`,
fuera de Git, antes de descargar. Ningún FASTQ crudo válido fue reemplazado.

## Validación previa al análisis

La auditoría del manifest registra 186/186 FASTQ válidos, 394.515.770 bytes,
cero pendientes y cero exclusiones nuevas. Se comprobaron tamaño, MD5,
gzip/CRC y estructura FASTQ completa. El informe y su procedencia están en
`results/planning/pacheco_montealegre_2020_colombia/download_status*.json`.
La ejecución de la auditoría registra `f908ba6` y conserva el estado Git observado.
La comprobación renovada de NCBI/ENA confirma los 94 runs sin discrepancias;
manifest y log de streams conservan su contenido. Las auditorías auxiliares
iniciales registran Python del sistema 3.14.4. Tras el ajuste operativo, el
checkpoint y la auditoría de los 186 archivos se repitieron sobre `0941180`
en el entorno Micromamba registrado, con Python 3.12.13 y las versiones fijadas.
El manifest y la planificación renovados pasan una comprobación independiente
de 31 entradas y seis artefactos, incluida la configuración operativa vigente.

Los 18 componentes del entorno existente cargan correctamente y conservan las
versiones registradas, incluidos R 4.4.3, DADA2 1.34.0, Cutadapt 5.0 y
Snakemake 8.30.0. Las cuatro regresiones de descarga pasan; el dry-run identifica
670 trabajos para completar las etapas pendientes. Los avisos de metadata
histórica ausente del scheduler no invalidan los FASTQ verificados.

## Alcance de la reanudación

Se conserva el inventario de 94 runs y el manifest paired-end de 93 runs,
cuatro lotes y 47 observaciones lote-tiempo. No cambian criterios de inclusión,
primers, filtros, taxonomía, tratamiento de ceros ni semillas.

La siguiente ejecución debe completar QC, DADA2, taxonomía, separación
bacteriana y diversidad descriptiva con validación de cada etapa. Este registro
de descarga no constituye un resultado biológico, un contraste temporal ni un
meta-análisis. El piloto de seis muestras sigue siendo la única ejecución
científica completa hasta que se cierren y validen esas etapas.

## Recursos operativos

La suite completa aprobó 144 pruebas en el entorno registrado. La primera
ejecución inició los 670 trabajos con cuatro tareas concurrentes y cerró las
validaciones de descarga y estructura de los 186 FASTQ.

Se verificaron ocho procesadores y 5.875 MB de RAM en WSL, con 4.845 MB
disponibles en ese momento. Se interrumpió ordenadamente el scheduler para
ampliar el límite a ocho tareas; el presupuesto global sigue en 4.500 MB y
DADA2/taxonomía siguen usando dos hilos. Ocho FastQC reservan 4.096 MB.
Los resultados FastQC completos se conservan y el reinicio revisa trabajos
incompletos. No se modifican entradas ni parámetros científicos.

## Composición descriptiva del piloto

La política `61da448` y la implementación `ad4c3ca` añaden resúmenes por
familia/género, conservando linajes y desconocidos. El perfil piloto produce
95 filas de composición para las seis muestras, 88 ASVs y 123.693 lecturas.
La validación reconcilia cada linaje y cobertura, comprueba 15 entradas y cinco
artefactos y se repite con `--validate-only` sin regenerar salidas. Las siete
pruebas nuevas usan datos sintéticos; los resultados anteriores proceden de
las muestras públicas reales. No se agregan estratos ni se contrastan hipótesis.
La suite completa posterior aprueba 151 pruebas en Linux (580,497 segundos
con el procesamiento del estudio en paralelo), incluidas las integraciones
existentes de R y Snakemake. No hay nuevas dependencias.
