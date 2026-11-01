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

## Diagnóstico inicial del estudio completo

Los 186 informes FastQC crudos cubren 2.728.920 pares, con 4.955–79.920
pares por muestra. La calidad por base pasa en 184 archivos y muestra avisos
en los dos mates de `SRR7899641` (`BactSaII120z2`). Los 186 archivos pasan
calidad por secuencia.

El módulo de adaptadores pasa en 183 archivos: ambos mates de `SRR7899641`
tienen estado `fail` y R2 de `SRR7899691` (`BactSaI096z1`) tiene `warn`.
Son estados de FastQC, no exclusiones ni una prueba de contaminación biológica.
La tabla fuente es `results/qc/pacheco_montealegre_2020_colombia/full/raw_read_quality.tsv`.
Se revisarán de nuevo tras el recorte y junto a la retención/longitud DADA2;
no se modifican los parámetros registrados para eliminar estas alertas.

La profundidad desigual debe considerarse antes de interpretar riqueza y efectos
temporales. Los 93 runs siguen representando cuatro lotes independientes y
47 observaciones lote-tiempo, con sus estratos separados.

## Corrección del cierre QC

El primer intento cerró 622 trabajos y se detuvo al resumir Cutadapt: el
validador compartido derivaba siempre una salida `pilot`, aunque el ejecutor
había escrito correctamente los 93 runs en el directorio `full` registrado.
No hubo un cambio de parámetros ni una pérdida de archivos crudos.

Se permite pasar el directorio de procesamiento registrado al validador;
permanece la ruta piloto por defecto. Se mantienen las comprobaciones de
entradas, estudio, run, mate, primers, flags y conteos. Se rechazan directorios
fuera del estudio y rutas de salida de otro run. Una nueva integración ejecuta
Cutadapt 5.0 con FASTQ sintéticos en la ruta completa, sin sustituir el parser
por un mock. El resumen de los 93 informes reales ya valida correctamente.
La reanudación renueva las etapas cuyo código de entrada cambió y conserva
los crudos protegidos.
La corrección se registra en `67c3820`; la suite completa pasa 154 pruebas
en 101,560 segundos antes de reanudar. No se alteran los parámetros científicos.

## Cierre QC del estudio completo

El QC renovado cierra sobre `1250c72`: 505 entradas y 11 artefactos pasan
la validación independiente. Los 93 runs conservan los 2.728.920 pares después
de Cutadapt (100 %), sin exclusiones. Se recortaron tres R1 y 109 R2, retirando
97 y 3.715 bases respectivamente. Son recortes de primers bajo la política
registrada; no se añadió recorte general de adaptadores o calidad.

Los 186 controles posteriores están completos: 184 pasan calidad por base y
los dos mates de SRR7899641 conservan `warn`. Adaptadores mantiene 183 `pass`,
dos `fail` (SRR7899641) y un `warn` (SRR7899691 R2). Estos avisos se mantienen
para revisar retención y secuencias DADA2; la validación de integridad y balances
no equivale a afirmar que todos los módulos FastQC pasan.

Los resúmenes, configuración y procedencia están en
`results/qc/pacheco_montealegre_2020_colombia/full/`; los informes completos y
FASTQ permanecen fuera de Git. El workflow inicia DADA2 con los mismos
parámetros científicos registrados. Las marcas de tiempo de los outputs se
conservan en UTC; este reporte utiliza la fecha local de Colombia.
