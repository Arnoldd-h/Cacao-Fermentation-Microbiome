# Scripts

Los entry points actuales están en `metadata/`, incluida la comparación
reproducible del dataset piloto. Los módulos futuros se añadirán
por responsabilidad (`download`, `qc`, `dada2`, `taxonomy`, `harmonization`,
`diversity`, `temporal`, `differential_abundance`, `meta_analysis`, `figures`)
cuando exista implementación real, evitando directorios vacíos y placeholders.

`metadata/download_pilot_fastq.py` es reanudable mediante archivos `.part`, no
sobrescribe un FASTQ final que falle validación y exige coincidencia de bytes y
MD5 antes de promover una descarga a archivo final.
