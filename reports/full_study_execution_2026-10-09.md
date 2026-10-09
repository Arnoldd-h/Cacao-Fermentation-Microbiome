# Avance de la ampliación — 2026-10-09

El workflow ampliado está implementado y probado. La descarga conserva todos
los archivos disponibles y el análisis completo espera un FASTQ público;
este reporte no presenta resultados DADA2, taxonómicos ni de diversidad del
estudio completo como si ya hubieran sido ejecutados.

| Etapa | Estado observado |
|---|---|
| NCBI/ENA | 94 runs incluidos verificados; ninguna discrepancia de metadata |
| Elegibilidad paired-end | 93 runs, cuatro lotes, 47 observaciones lote-tiempo |
| Streams públicos | 189 registrados, 186 seleccionados, tres aislados excluidos |
| Descarga local | 185/186 FASTQ válidos; 393.057.692 bytes |
| Integridad | Bytes, MD5, gzip/CRC y estructura FASTQ comprobados íntegramente |
| QC/DADA2 completos | Pendientes: falta R2 de SRR7899804 |
| Pruebas | 144 pruebas Python aprobadas, incluidas 19 de ampliación |
| Piloto previo | 20 entradas/17 artefactos válidos; siete TSV sin cambios por SHA-256 |

`SRR7899803` conserva su inclusión en el inventario, pero no tiene pares
completos y quedó excluido de la ruta paired-end antes de analizar resultados.
`SRR7899804` sí permanece en el manifest: la API declara un R2 de 1.458.078
bytes y MD5 `49b5d1b26dd75169352b52aace2c5524`, pero el enlace devuelve un
directorio HTML vacío. Son dos decisiones distintas; no se excluye una muestra
por un fallo de transferencia ni se relaja la validación para continuar.

Las comprobaciones HEAD abarcaron los 186 enlaces; sólo el R2 indicado
devolvió HTML. Se verificaron HTTP/HTTPS, un segundo hostname público y FTP,
además del registro ENA. Dos reintentos adicionales de la ejecución no
recuperaron ese archivo. Una transferencia detenida de SRR7899805 se reanudó
con Range y su archivo final pasó MD5; no se alteraron archivos crudos válidos.

Los registros pequeños están en
`results/planning/pacheco_montealegre_2020_colombia/`: `source_audit.json`,
`fastq_stream_log.tsv`, `resource_plan.json`, `packaging_source_audit.json`,
`packaging_audit.json` y `download_status.json`, con snapshots/hashes.
El manifest procede de `fa236ed`; la auditoría de descarga, de `cf99de0`.
Los estados Git registrados incluyen los nuevos resultados aún sin commit.
La renovación del piloto procede de `e472984` y reproduce sus siete TSV.

No se publica un tag de análisis completo ni el hito multistudio. La próxima
acción es recuperar ese R2 desde el origen público o registrar y validar una
ruta alternativa que preserve secuencias, calidades e identificadores, antes
de ejecutar el estudio de 93 runs. Después se reanuda el mismo workflow;
los FASTQ verificados y SILVA protegida se reutilizan.
