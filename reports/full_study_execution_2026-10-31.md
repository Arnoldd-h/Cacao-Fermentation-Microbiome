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

## Cierre DADA2 del estudio completo

La ejecución sobre `1250c72` conserva las 93 muestras: 2.540.412 pares filtrados,
2.395.174 ensamblados y 2.254.530 sin quimeras (82,6162 % de 2.728.920 de
entrada). De 1.606 secuencias ensambladas se conservan 782 ASVs después de
identificar 824 quimeras por consenso. Todavía no son una tabla bacteriana
separada por taxonomía. Las longitudes retenidas van de 236 a 247 nt;
no se aplica un filtro de longitud no registrado.

Los modelos F/R convergen en seis rondas, antes del límite de diez. Los
diagnósticos registran 101.383.928 y 101.237.542 observaciones de transición.
Se revisaron ambos PNG: curvas y puntos legibles; se conservan nueve warnings
de representación logarítmica de frecuencias cero, sin ocultarlos. La
convergencia no equivale a ausencia de errores o contaminación.

El validador independiente comprueba 192 entradas y 18 artefactos, conteos,
secuencias, longitudes y modelos. Los hashes de entrada se verificaron también
al cerrar R. La procedencia conserva `git_dirty=true` por resultados QC y
planificación aún sin commit al iniciar, cuyos contenidos se registran por hash.
Los cambios posteriores de documentación y productores descriptivos no alteraron
las entradas científicas durante DADA2.

SRR7899641 conserva 9.178 pares sin quimeras de 22.820 (40,22 %); SRR7899691
conserva 28.413 de 33.607 (84,54 %). Los avisos FastQC permanecen documentados,
sin excluir estas muestras ni modificar parámetros. Se inicia la clasificación
SILVA con el mismo umbral primario y sensibilidad registrados para el piloto.

## Revisión de figuras y renovación de procedencia

El primer recorrido completo cerró sus 481 trabajos y validó taxonomía,
separación y diversidad. La revisión visual de cobertura detectó un título
fijo `Pilot coverage` en el productor compartido. Se sustituye por
`Taxonomic coverage`, válido para ambos alcances, sin modificar parámetros
ni cálculos. Se guardan hashes de las tablas anteriores para contrastar
que la renovación conserve sus valores. Las figuras se regeneran por código;
no se retocan manualmente.

La renovación sobre `da43eeb` cierra siete trabajos y vuelve a validar todas
las etapas. La comparación encuentra idénticas la taxonomía primaria y las
once tablas bacterianas/de diversidad. Cambian soportes y sensibilidad a 50,
con una diferencia neta de una ASV/cinco lecturas asignadas a género. El
desempate nativo no está completamente controlado por la semilla R; evidencia,
fuentes y limitación pendiente en `taxonomy_reproducibility_2026-10-31.md`.

## Cierre taxonómico del estudio completo

SILVA 138.2 conserva las 782 ASVs y 2.254.530 lecturas originales; todas tienen
Kingdom asignado como Bacteria. La máscara primaria a 80 asigna género a
475 ASVs y 1.855.013 lecturas, y familia a 692 ASVs y 1.932.980 lecturas.
La última sensibilidad a 50 asigna género a 648 ASVs y 2.051.009 lecturas.
Se marcan 18 ASVs con etiquetas de orgánulos para la separación posterior;
no se eliminan de estas tablas originales y no hay asignación de especie.

La validación independiente comprueba 21 entradas y 15 artefactos, sin
warnings R en esta clasificación. El PNG de cobertura se revisó con su título
corregido y ambas sensibilidades; se conservan también PDF y SVG. Ausencia
de warnings R no elimina la limitación del RNG nativo documentada aparte.

## Separación y diversidad completas

La política primaria separa 18 ASVs etiquetadas como orgánulos y 157.550
lecturas, conservando 764 ASVs, 2.096.980 lecturas y todas las 93 muestras.
La separación pasa 21 entradas/nueve artefactos. No se añadieron exclusiones
de muestras, filtros de abundancia o umbrales de profundidad.

La diversidad descriptiva pasa 20 entradas/17 artefactos: riqueza de 1–75 ASVs,
Shannon efectivo de 1–30,363 y profundidad de 986–69.751 lecturas. Se conservan
CLR/Aitchison con pseudoconteos 1 y 0,5, Bray-Curtis y PCA; no hay rarefacción
ni pruebas inferenciales. Las primeras dos componentes primarias representan
el 12,68 % de la variación CLR. No hay warnings R de diversidad.

Las figuras alfa y PCA se revisaron visualmente: alcance completo, leyendas,
ejes, sensibilidad y ausencia de 93 etiquetas superpuestas. Alfa usa índice
de muestra, no horas; las figuras temporales posteriores usan las horas reales.
Los once TSV bacterianos/de diversidad son idénticos antes y después de
renovar el título taxonómico. El CLR largo de 16,3 MB queda ignorado y se
regenera mediante el workflow; su hash y validación se conservan.

La revisión de bytes detectó finales de línea mixtos en la copia local de
`config/full_study_diversity.yaml`, aunque el blob Git registrado ya usa LF.
Se normalizó la copia local a LF y se comprobó igualdad JSON de todos sus
valores, además de coincidencia exacta con el blob Git
`8e761c3f4b5eb528c7e22beb5437bc4f1e08d863`. No hay un cambio metodológico ni
un diff de configuración. Se renovaron cinco trabajos de separación/diversidad
y sus resúmenes para que snapshots y hashes correspondan a esos mismos bytes
en un checkout nuevo. Se mantienen idénticos los once TSV derivados y cuatro
productos de composición/representación, incluido el PNG temporal revisado.

## Renovación del piloto y cierre de workflows

El workflow piloto cierra 33 trabajos después de los cambios compartidos de
QC y del título taxonómico. Sus ocho tablas DADA2, ocho taxonómicas y once
bacterianas/de diversidad son idénticas a la referencia anterior. Los cinco
registros QC pasan 257 comprobaciones independientes de bytes y SHA-256.
El piloto conserva 93 ASVs/156.561 pares sin quimeras y 88 ASVs/123.693
lecturas después de separar orgánulos. Se renuevan también composición y
figuras temporales; el PNG temporal es idéntico al ya revisado.

Los dry-runs finales de `workflow/Snakefile` y `workflow/study.Snakefile`
terminan con `Nothing to be done`. Snakemake conserva avisos sobre metadata
auxiliar histórica ausente en 12 descargas piloto y 20 trabajos completos
de descarga/FastQC. No son trabajos pendientes: sus archivos y procedencia
científica siguen comprobados, sin reemplazar crudos para reconstruir esa caché.

## Composición, figuras temporales y alcance del hito

El perfil completo contiene 1.477 filas de composición por familia/género,
93 muestras, 764 ASVs y todas las 2.096.980 lecturas retenidas. Pasa 15
entradas/cinco artefactos al generar y con `--validate-only`. Se conservan
linajes, etiquetas vacías y denominadores completos por muestra.

Las figuras temporales describen ocho series de lote/estrato de los cuatro
lotes, con horas originales, riqueza, Shannon efectivo y profundidad. Pasan
16 entradas/ocho artefactos al generar y con `--validate-only`; PNG revisado
visualmente y exportaciones PDF/SVG/PNG de 300 dpi. No se imputan huecos ni
se conectan lotes o estratos distintos.

Las notas de resultados registran la fracción etiquetada como Acetobacter a
0 horas (0,14–14,42 %) y en el último tiempo observado (58,27–99,58 %), para
las ocho series. Son observaciones composicionales, con cobertura y profundidad
desiguales, cuatro lotes independientes y últimos tiempos distintos (108–144 h).
No constituyen pruebas longitudinales, crecimiento absoluto, core ni evidencia
de conservación entre estudios. El hito descriptivo conserva la limitación
bootstrap y las tareas inferenciales pendientes.
