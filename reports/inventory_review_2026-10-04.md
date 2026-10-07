# Revisión del inventario iniciada el 4 de octubre de 2026

Revisión iniciada el 2026-10-04 UTC y completada el 2026-10-07 UTC
(6 de octubre en Colombia). Las fechas reales de recuperación quedan en los
TSV. No se descargan secuencias en esta revisión.

## Evidencia y decisiones consolidadas

- La consulta NCBI configurada recuperó 34 BioProjects. Se revisaron registros
  ENA de los 34, publicaciones primarias, XML de experimentos y nombres de los
  archivos depositados; estos nombres son metadata, no lecturas descargadas.
- PRJNA865318: la publicación primaria DOI
  [10.1094/PHYTOFR-08-23-0104-R](https://apsjournals.apsnet.org/doi/abs/10.1094/PHYTOFR-08-23-0104-R)
  y el [repositorio de autores](https://github.com/nnguyenlab/cacao-microbiome)
  resuelven el pendiente: las corridas bacterianas de fermentación depositadas
  sólo representan el día 1. Las filas planificadas de días 3/6 en la hoja de
  autores no prueban la existencia de FASTQ de esos tiempos. Se excluye de la
  serie longitudinal primaria.
- PRJNA1104253: [métodos primarios](https://pmc.ncbi.nlm.nih.gov/articles/PMC12408344/)
  confirman que el componente natural on-farm es WGS; los amplicones proceden
  de ensayos controlados con manipulación aséptica, temperatura programada y,
  para dropout, granos esterilizados/pulpa artificial. Los primers 338F/806R
  quedan verificados. Se conserva separado del análisis primario natural.
- PRJEB40850: [artículo primario](https://doi.org/10.3389/fmicb.2020.616875)
  identifica F1/F2 como controles espontáneos de 20 kg, noviembre de 2017,
  92 horas. Los nombres ENA `submitted_ftp` terminados en `_16S.fastq.gz`
  permiten aislar 12 runs V4 frente a ITS. No se infiere el marcador desde el
  tamaño del archivo.
- PRJEB57747: [artículo primario](https://doi.org/10.3389/fmicb.2023.1232323)
  identifica F01/F02 como dos fermentaciones espontáneas de 36 kg, noviembre
  de 2019, 120 horas. Sus 16 runs Sequel II se depositaron con sufijo
  `.16S.fastq.gz`. El par 27F/1492R se verifica en el método explícitamente
  citado, [Callahan et al. 2019](https://doi.org/10.1093/nar/gkz569).
- PRJNA420946: el artículo [Mota-Gutierrez et al. 2018](https://pmc.ncbi.nlm.nih.gov/articles/PMC6147003/)
  enlaza explícitamente SRP126069, cuyo registro ENA pertenece a PRJNA420946.
  Declara fin de fermentación a 120 h, mientras ENA contiene controles a
  144 h. No se eliminan ni recodifican esas observaciones para forzar
  concordancia: el candidato permanece pendiente hasta aclarar el conflicto.
- PRJEB53853: el artículo [Díaz-Muñoz et al. 2023](https://doi.org/10.1016/j.fm.2022.104115)
  resuelve los tratamientos F1/F2, pero las bibliotecas depositadas tienen
  descripción conjunta 16S/ITS y nombres enviados sin etiqueta de marcador.
  Hace falta demostrar una separación de lecturas por marcador; el artículo
  informa fracciones bacterianas y fúngicas dentro de 60 muestras. Se mantiene
  pendiente, sin presentar todos los reads como 16S.

## Resultado verificable y procedencia

Los 34 BioProjects tienen decisión documentada: 4 incluidos, 7 pendientes y
23 excluidos. Se recuperaron 2.357 runs: 182 incluidos, 93 pendientes y 2.082
excluidos. `samples.tsv` contiene los 275 runs candidatos incluidos o
pendientes. El log tiene 2.112 decisiones: 2.105 exclusiones (2.082 runs y
23 estudios) y 7 pendientes a nivel estudio; no son 2.112 muestras excluidas.

El conjunto incluido reúne 12 etiquetas de lotes: Colombia 4, México 4,
Costa Rica 2017 2 y Costa Rica 2019 2. Las 182 corridas no son 182 unidades
experimentales independientes. Sus FASTQ suman 2.056.243.128 bytes estimados
por ENA. La suma de todos los proyectos, incluyendo WGS y otras exclusiones,
es 702.333.578.305 bytes; son estimaciones, no archivos descargados.

Se aplicó la skill local `sra-metadata-inventory`. Se consultaron NCBI BioProject,
la API ENA para todos los runs y XML de experimentos cuando el marcador no
estaba resuelto. Se buscaron accessions, títulos exactos y combinaciones de
autores/variedades en publicaciones primarias, Europe PMC y repositorios de
autores. Los textos completos accesibles de PMC, Frontiers y UGent permitieron
revisar métodos y referencias a protocolos. Los cachés locales quedan fuera de
Git; las fuentes duraderas están en este informe y `config/datasets.yaml`.

El primer descubrimiento de la sesión final recibió HTTP 500 de ENA para
PRJEB101075 tras tres intentos y falló sin publicar una tabla parcial. Una
nueva ejecución completa terminó con los 34 registros. También terminó la
reconstrucción completa y pasó el validador independiente. La fecha exacta
de cada recuperación está en `retrieved_at_utc`. El commit base al validar fue
`dd4f3649a823afef9d0a07953c05711625c64d04`, con cambios aún sin commit: ese hash
no identifica por sí solo el inventario revisado.

`submitted_ftp` y `submitted_format` conservan metadata observada de los archivos
enviados. Sus nombres distinguen marcadores en Costa Rica 2017/2019, pero no
prueban ausencia de procesamiento previo de reads. El layout procede del
registro de biblioteca. No se usa tamaño de archivo como prueba de marcador.
Los valores desconocidos quedan vacíos; las descripciones anteriores
`unknown`/`not reported` permanecen en notas de configuración. Un estudio no
sustituye un lote desconocido y el año del depósito no se toma como publicación.

## Detalle de los nuevos incluidos

**PRJEB40850, Costa Rica 2017.** Los métodos y la figura del
[artículo primario](https://pmc.ncbi.nlm.nih.gov/articles/PMC7829357/)
identifican F1/F2 como controles espontáneos de 20 kg de cacao Trinitario,
fermentados en recipientes plásticos en Heredia en noviembre de 2017. Los seis
tiempos son 0/6/20/44/68/92 h. El selector exige AMPLICON, PAIRED, F1/F2 y un
archivo enviado `_16S.fastq.gz`. El protocolo explícitamente citado,
[De Bruyn et al. 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5165123/),
reporta los segmentos no-plataforma F515 `GTGTGCCAGCMGCCGCGGTAA` y R806
`GGACTACHVGGGTWTCTAAT`. Se conserva el `GT` inicial reportado y se anotan las
construcciones completas. Su presencia efectiva en los reads deberá verificarse
en un piloto propio antes de parametrizar este estudio.

**PRJEB57747, Costa Rica 2019.** Los
[métodos primarios](https://pmc.ncbi.nlm.nih.gov/articles/PMC10445768/)
identifican dos controles espontáneos F01/F02 de 36 kg, fermentados en Heredia
en noviembre de 2019. Sus ocho tiempos son 0/6/12/24/48/72/96/120 h. El selector
exige AMPLICON, SINGLE, Sequel II, F01/F02 y `.16S.fastq.gz`, separando ITS MiSeq
y WGS NovaSeq. El método citado,
[Callahan et al. 2019](https://pmc.ncbi.nlm.nih.gov/articles/PMC6765137/),
verifica 27F `AGRGTTYGATYMTGGCTCAG` y 1492R `RGYTACCTTGTTACGACTT`.
La aceptación científica no habilita automáticamente el workflow Illumina para
PacBio: se requiere un procesamiento específico.

Los experimentos de 2017/2019 comparten equipo, sitio y procedencia del cacao,
pero se realizaron en años distintos. Se mantienen estudios y lotes separados;
la integración deberá evaluar dependencia por laboratorio/sitio. El criterio
vigente no exige literalmente fermentar en una finca: los controles espontáneos
de masa real cumplen los criterios documentados, conservando ese contexto.

## Los siete pendientes y el requisito para resolverlos

| BioProject | Runs candidatos | Requisito pendiente |
|---|---:|---|
| PRJEB53853 | 0 | Separación validada de los marcadores en bibliotecas mixtas 16S/ITS |
| PRJEB82327 | 48 | Protocolo primario vinculado al acceso, primers e independencia de cohorte |
| PRJEB82871 | 16 | Protocolo primario vinculado al acceso, primers e independencia de cohorte |
| PRJNA420946 | 14 | Reconciliar 120 h del artículo con 144 h depositadas y lotes incompletos |
| PRJNA407677 | 4 | Mapa run-tiempo-lote, protocolo primario y primers |
| PRJNA935329 | 6 | Protocolo original, diseño espontáneo, lote, duración y primers |
| PRJNA842267 | 5 | Métodos completos, manejo del material y unidades del control F0 |

**Costa Rica 2018:** el texto completo se recuperó del
[repositorio UGent](https://backoffice.biblio.ugent.be/download/01HWCNPM6MZCW83NM90MXTM5VZ/01HWCNRFZKZ1B7FMS24E3P0X04).
F1/F2 y 96 h quedaron resueltos, pero los nombres enviados no separan 16S/ITS.
No se presentan las 60 bibliotecas como bacterianas. Las secuencias exactas de
primers quedan vacías hasta verificar el protocolo completo citado y la
separación de marcadores.

**Ecuador Nacional/CCN-51:** ENA identifica 48/16 runs full-length 16S,
respectivamente, con ocho tiempos y dos lotes por proyecto. Las búsquedas por
accession, título y autores no recuperaron artículo/protocolo primario con
primers y vínculo inequívoco entre acceso y método. El
[resumen doctoral de Dario Van de Voorde](https://www.vub.be/sites/default/files/2025-05/PhD%20template%20ENG%20Dario%20Van%20de%20Voorde.pdf)
describe experimentos concordantes, pero no acredita primers ni independencia
de cohortes. Se necesita la tesis completa o artículo con accessions/métodos.

**Camerún:** [SRX3441502](https://www.ebi.ac.uk/ena/browser/api/xml/SRX3441502)
verifica la relación SRP126069–PRJNA420946. SRR6342740/SRR6342743 declaran 144 h,
frente al final de 120 h del artículo. Se retienen esos tiempos y `pending`;
duración, tiempo relativo y etapa quedan vacíos. Hay tres etiquetas de lotes
de control: heap primero/segundo y box segundo; el depósito no reconstruye
todos los duplicados descritos. Los primers V3–V4 se verificaron mediante
[Klindworth et al. 2013](https://pmc.ncbi.nlm.nih.gov/articles/PMC3592464/),
citado en métodos. Se necesita una hoja original o corrección del depósito que
explique 144 h y las unidades faltantes; no se eliminan observaciones para
forzar concordancia.

**Brasil PRJNA407677:** los XML de
[SRX3270869](https://www.ebi.ac.uk/ena/browser/api/xml/SRX3270869),
[SRX3274478](https://www.ebi.ac.uk/ena/browser/api/xml/SRX3274478),
[SRX3274479](https://www.ebi.ac.uk/ena/browser/api/xml/SRX3274479) y
[SRX3274480](https://www.ebi.ac.uk/ena/browser/api/xml/SRX3274480)
identifican Ilheus/Tomeacu/Medicilandia/Placas 16S frente a dos ITS. La
descripción común menciona días 2/4/6 sin asignarlos a runs y no demuestra si
hubo mezcla de tiempos. El artículo
[Serra et al., DOI 10.1016/j.lwt.2019.02.038](https://doi.org/10.1016/j.lwt.2019.02.038)
es una pista, no un vínculo verificado con el accession: el editor negó el
texto completo con HTTP 403. No se rellenan tiempos, lotes o primers mediante
una reanálisis secundaria.

**Tabasco PRJNA935329:** el XML
[SRX19376241](https://www.ebi.ac.uk/ena/browser/api/xml/SRX19376241)
confirma V3–V4 bacteriano. Day1–Day6 permiten conservar 24–144 h, pero no
prueban lote ni duración. Las búsquedas recuperan una reanálisis geográfica
posterior, no el protocolo experimental original. Los campos no verificados
quedan vacíos; seis días no se convierten en seis fermentaciones independientes.

**Cocktails PRJNA842267:** el resumen primario de
[Investigating luxS gene expression in lactobacilli along lab-scale cocoa fermentations](https://doi.org/10.1016/j.fm.2023.104429)
distingue F0 no inoculado de F1–F3 con cocktails. El XML
[SRX15471927](https://www.ebi.ac.uk/ena/browser/api/xml/SRX15471927)
proporciona primers V3–V4. F0 tiene 0/24/48/72/96 h sin etiqueta extra de réplica,
a diferencia de los tratamientos. Se recuperó el resumen mediante Europe PMC;
el editor y el [preprint bioRxiv](https://doi.org/10.1101/2022.06.14.496151)
no permitieron recuperar métodos completos. Falta verificar manejo del material,
comparabilidad del control y correspondencia con unidades independientes. Los
cinco tiempos F0 no equivalen a cinco réplicas biológicas.

## Siete exclusiones de la cola manual

- **PRJEB101075:** el proyecto primario declara ITS full-length y comunidades
  fúngicas mock de cacao; los alimentos naturales son sourdough/lambic.
- **PRJEB50472:** el [artículo primario RSC](https://pubs.rsc.org/en/content/articlehtml/2021/fo/d1fo01155c)
  estudia digestión-fermentación gastrointestinal in vitro, no masa pulpa-grano.
  Se conservan las descripciones depositadas de suelo, que son discordantes;
  no se corrigen inventando etiquetas.
- **PRJNA420973:** el diseño de
  [SRX3442157](https://www.ebi.ac.uk/ena/browser/api/xml/SRX3442157)
  verifica ITS2, compañero fúngico de Camerún.
- **PRJNA475867:** el diseño de
  [SRX4201987](https://www.ebi.ac.uk/ena/browser/api/xml/SRX4201987)
  declara 26S de levaduras, aunque la estrategia agregada sea AMPLICON.
- **PRJNA783055:** las 52 descripciones de run identifican contenido cecal
  murino tras suplementación dietaria con cacao; la matriz es intestino.
- **PRJNA842340:** el proyecto declara el componente ITS de PRJNA842267;
  títulos y alias son concordantes.
- **PRJNA935525:** el diseño de
  [SRX19386157](https://www.ebi.ac.uk/ena/browser/api/xml/SRX19386157)
  verifica ITS, compañero fúngico de PRJNA935329.

## Tabla completa de decisiones

«Candidatos» cuenta los runs incluidos o pendientes que cumplen el selector.
Cero en un pendiente puede indicar falta de selector seguro, no ausencia de
reads. Cada accession enlaza la fuente primaria ENA; los artículos específicos
figuran arriba y en configuración.

| BioProject / registro primario | Decisión | Runs | Candidatos | Evidencia decisiva |
|---|---|---:|---:|---|
| [PRJNA492720](https://www.ebi.ac.uk/ena/browser/view/PRJNA492720) | incluir | 292 | 94 | Serie espontánea V4 verificada |
| [PRJNA865318](https://www.ebi.ac.uk/ena/browser/view/PRJNA865318) | excluir | 63 | 0 | Sólo día 1 bacteriano recuperado |
| [PRJNA1104253](https://www.ebi.ac.uk/ena/browser/view/PRJNA1104253) | excluir | 473 | 0 | Amplicones de ensayos controlados; natural es WGS |
| [PRJNA552479](https://www.ebi.ac.uk/ena/browser/view/PRJNA552479) | excluir | 14 | 0 | WGS |
| [PRJNA1257864](https://www.ebi.ac.uk/ena/browser/view/PRJNA1257864) | excluir | 10 | 0 | WGS |
| [PRJNA1264670](https://www.ebi.ac.uk/ena/browser/view/PRJNA1264670) | excluir | 9 | 0 | WGS |
| [PRJNA627078](https://www.ebi.ac.uk/ena/browser/view/PRJNA627078) | incluir | 132 | 60 | Serie espontánea V3–V4; cuatro cajas |
| [PRJEB53853](https://www.ebi.ac.uk/ena/browser/view/PRJEB53853) | pendiente | 60 | 0 | Bibliotecas 16S/ITS sin separación validada |
| [PRJEB82327](https://www.ebi.ac.uk/ena/browser/view/PRJEB82327) | pendiente | 369 | 48 | Protocolo/primers no verificados |
| [PRJEB82871](https://www.ebi.ac.uk/ena/browser/view/PRJEB82871) | pendiente | 140 | 16 | Protocolo/primers no verificados |
| [PRJNA962540](https://www.ebi.ac.uk/ena/browser/view/PRJNA962540) | excluir | 11 | 0 | Sin serie 16S temporal verificable |
| [PRJEB40850](https://www.ebi.ac.uk/ena/browser/view/PRJEB40850) | incluir | 108 | 12 | F1/F2 espontáneos, 92 h, V4 |
| [PRJEB57747](https://www.ebi.ac.uk/ena/browser/view/PRJEB57747) | incluir | 243 | 16 | F01/F02 espontáneos, 120 h, full-length |
| [PRJNA420946](https://www.ebi.ac.uk/ena/browser/view/PRJNA420946) | pendiente | 42 | 14 | Conflicto 120/144 h y lotes incompletos |
| [PRJNA407677](https://www.ebi.ac.uk/ena/browser/view/PRJNA407677) | pendiente | 6 | 4 | Tiempo/lote y primers no verificables por run |
| [PRJNA935329](https://www.ebi.ac.uk/ena/browser/view/PRJNA935329) | pendiente | 6 | 6 | Protocolo/primers/lote pendientes |
| [PRJNA842267](https://www.ebi.ac.uk/ena/browser/view/PRJNA842267) | pendiente | 35 | 5 | Métodos completos del control F0 pendientes |
| [PRJEB101075](https://www.ebi.ac.uk/ena/browser/view/PRJEB101075) | excluir | 25 | 0 | ITS fúngico; cacao mock |
| [PRJEB50472](https://www.ebi.ac.uk/ena/browser/view/PRJEB50472) | excluir | 21 | 0 | Modelo gastrointestinal |
| [PRJNA420973](https://www.ebi.ac.uk/ena/browser/view/PRJNA420973) | excluir | 41 | 0 | ITS2 |
| [PRJNA475867](https://www.ebi.ac.uk/ena/browser/view/PRJNA475867) | excluir | 28 | 0 | 26S levaduras |
| [PRJNA783055](https://www.ebi.ac.uk/ena/browser/view/PRJNA783055) | excluir | 52 | 0 | Contenido cecal murino |
| [PRJNA842340](https://www.ebi.ac.uk/ena/browser/view/PRJNA842340) | excluir | 35 | 0 | ITS |
| [PRJNA935525](https://www.ebi.ac.uk/ena/browser/view/PRJNA935525) | excluir | 6 | 0 | ITS |
| [PRJDB13568](https://www.ebi.ac.uk/ena/browser/view/PRJDB13568) | excluir | 12 | 0 | WGS |
| [PRJEB38017](https://www.ebi.ac.uk/ena/browser/view/PRJEB38017) | excluir | 16 | 0 | RNA-Seq/WGS |
| [PRJEB46227](https://www.ebi.ac.uk/ena/browser/view/PRJEB46227) | excluir | 2 | 0 | WGS |
| [PRJEB57751](https://www.ebi.ac.uk/ena/browser/view/PRJEB57751) | excluir | 8 | 0 | WGS |
| [PRJNA1002093](https://www.ebi.ac.uk/ena/browser/view/PRJNA1002093) | excluir | 21 | 0 | RNA-Seq |
| [PRJNA1037292](https://www.ebi.ac.uk/ena/browser/view/PRJNA1037292) | excluir | 21 | 0 | WGS |
| [PRJNA212389](https://www.ebi.ac.uk/ena/browser/view/PRJNA212389) | excluir | 24 | 0 | RNA-Seq |
| [PRJNA527768](https://www.ebi.ac.uk/ena/browser/view/PRJNA527768) | excluir | 10 | 0 | WGS |
| [PRJNA627182](https://www.ebi.ac.uk/ena/browser/view/PRJNA627182) | excluir | 20 | 0 | WGS |
| [PRJNA83439](https://www.ebi.ac.uk/ena/browser/view/PRJNA83439) | excluir | 2 | 0 | WGS |

Las diez exclusiones anteriormente automáticas se revisaron contra todos los
runs ENA y se configuraron explícitamente: siete WGS, dos RNA-Seq y una mezcla
RNA-Seq/WGS. Ninguna contiene un subconjunto AMPLICON 16S. La decisión se apoya
en estrategia depositada, no sólo en título.

## Diseño, limitaciones y validación

La unidad independiente es la fermentación/lote. En México, tres extracciones
por caja-tiempo son submuestras; en Colombia, los estratos pertenecen al mismo
evento. Deben conservarse y agregarse o modelarse conforme al protocolo antes
de inferencia. El piloto de seis runs valida procesamiento, no efectos
multicéntricos ni réplicas biológicas por tiempo.

No se concatenarán ASVs V4, V3–V4 y full-length entre estudios. La integración
requiere taxonomía armonizada, normalmente género/familia, conservando rangos
originales. País, laboratorio, año, plataforma y manejo pueden estar
confundidos. Los siete pendientes quedan fuera del primario hasta resolver los
requisitos descritos. No se contactó a autores; los intentos públicos fallidos
no se presentan como evidencia primaria verificada.

Se completaron `discover_candidates.py`, `build_inventory.py` y
`validate_metadata.py`. La suite Python final ejecutó 85 tests: 81 pasaron y
4 comprobaciones de integración Snakemake se omitieron en Windows por requerir
el entorno Linux declarado. Incluye separación 16S/ITS, campos enviados,
ausencia de lotes inventados, conservación de 144 h pendientes y escritura TSV
que conserva el mtime si los bytes son idénticos.

La validación integrada posterior ejecutada por el proceso principal aprobó
89 tests en Linux, incluidas las cuatro integraciones Snakemake y las pruebas
independientes de artefactos DADA2. Los fixtures sintéticos sólo comprueban
contratos de software y no se presentan como evidencia de fermentación.

El ranking recalculado en memoria conserva PRJNA492720 primero. El manifest
resultante es idéntico al actual: seis runs y 25.154.687 bytes. Esta revisión
no modificó el manifest. `paired_end` en `studies.tsv` resume el BioProject
entero y puede ser `mixed` en PRJEB57747, aunque sus 16 candidatos sean todos
`SINGLE` PacBio; la metadata por run preserva esa diferencia.
