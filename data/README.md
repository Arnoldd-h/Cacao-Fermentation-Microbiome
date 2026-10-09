# Datos

- `raw/`: FASTQ y fuentes inmutables; no se versionan.
- `interim/`: artefactos regenerables entre pasos; no se versionan.
- `processed/`: productos analíticos finales pequeños; se versionan sólo cuando
  aportan interpretación y conservan procedencia.

Los directorios se crean por el workflow cuando son necesarios.

`processed/pilot/bacterial/` contiene los conteos enteros, secuencias ASV,
taxonomía primaria y metadata del piloto tras la separación registrada en
`config/diversity.yaml`: 88 ASVs, 123.693 lecturas y seis muestras.
Los conteos retenidos, secuencias, etiquetas y metadata coinciden con sus
entradas originales; no se les añade el pseudoconteo de los análisis CLR.

Las 93 ASVs y 156.561 lecturas originales permanecen en
`results/dada2/pilot/` y `results/taxonomy/pilot/`. Las decisiones de las 93 ASVs
y los balances por muestra están en `results/filtering/pilot/`, con hashes
y validación independiente. La derivación se regenera con `pilot_bacterial_table`
o como dependencia de `pilot_diversity`; ninguna tabla se edita manualmente.
