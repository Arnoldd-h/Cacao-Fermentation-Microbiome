# Reglas modulares

`pilot_qc.smk` implementa FastQC y MultiQC crudos para cada FASTQ declarado en
el manifest piloto, más el resumen machine-readable validado. El inventario
permanece en `workflow/Snakefile`; las etapas posteriores se incorporarán como
reglas pequeñas por estudio.
