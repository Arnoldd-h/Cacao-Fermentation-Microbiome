# Datos

- `raw/`: FASTQ y fuentes inmutables; no se versionan.
- `interim/`: artefactos regenerables entre pasos; no se versionan.
- `processed/`: productos analíticos finales pequeños; se versionan sólo cuando
  aportan interpretación y conservan procedencia.

Los directorios se crean por el workflow cuando son necesarios.
