# Métodos de clasificación taxonómica del piloto

## Alcance y referencia

La etapa consume las 93 ASVs no quiméricas del piloto PRJNA492720, conservando
su `study_id` y sus secuencias originales. El piloto comprende seis muestras;
las tablas describen resolución taxonómica técnica y no prueban sucesión.

Se utiliza [SILVA NR99 138.2 formateada para DADA2](https://doi.org/10.5281/zenodo.14169026),
archivo `silva_nr99_v138.2_toGenus_trainset.fa.gz`, enlazado por los
[mantenedores de DADA2](https://benjjneb.github.io/dada2/training.html).
El archivo comprimido tiene 139.996.892 bytes y MD5
`1764e2a36b4500ccb1c7d5261948a414`, según el depósito. Su SHA-256 verificado es
`e0af8246b900a922ea3ea68e7cb1a5e70ceed4db84246b6ea9865eb132452093`.
La inspección completa de gzip/FASTA cuenta 452.055 secuencias y 646.360.937 bases.
Se guarda la base bajo `references/databases/`, ignorada por Git; el pequeño
registro `references/silva_138.2.provenance.json` conserva su procedencia.

La referencia fue formateada con DADA2 1.35.4 según el proveedor. La ejecución
usa el clasificador ya fijado DADA2 1.34.0; no se actualiza el entorno ni la
nomenclatura SILVA. Deben citarse tanto el depósito formateado como el proyecto
SILVA: Quast et al. (2013), [doi:10.1093/nar/gks1219](https://doi.org/10.1093/nar/gks1219).

## Clasificación y soporte

Los parámetros de `config/taxonomy.yaml` y la decisión científica se fijaron
en `319f9ad`, antes de observar asignaciones. Se aplica `dada2::assignTaxonomy`
con comprobación de complemento inverso, seis rangos originales de Kingdom
a Genus y semilla 20260819. El generador usa Mersenne-Twister, Inversion y
Rejection; el paralelismo efectivo se registra en la procedencia.

El clasificador exporta la mejor llamada sin máscara (`minBoot=0`) y los
soportes bootstrap. La llamada primaria exige 80 en todos los ancestros;
la sensibilidad descriptiva exige 50 sobre el mismo ajuste. No se reasigna
un rango inferior si un ancestro falta o no supera el umbral. La
[guía oficial](https://benjjneb.github.io/dada2/assign.html) describe ambos
umbrales; 50 permite revisar la pérdida de resolución de estas ASVs cortas,
sin escoger parámetros según los taxones esperados o resultados temporales.

Los nombres originales se preservan, incluidos `uncultured` y otros grupos
sin nombre de género. La bandera `named_genus` permite distinguirlos de una
etiqueta de género resuelta; no los reemplaza. No se realiza asignación de
especies con la región V4. No se concatenan ASVs de otros estudios.

## Revisión de secuencias y conservación

Se marcan etiquetas exactas Chloroplast/Mitochondria, Kingdom conocido
distinto de Bacteria y Kingdom sin asignación. Las marcas se basan en la
llamada primaria y no excluyen ASVs ni lecturas. Una falta de asignación no
demuestra ausencia de eucariotas: esta referencia está optimizada para
Bacteria/Archaea. No hay controles negativos verificados y no se hace aquí
una evaluación de contaminación o abundancia diferencial.

La cobertura se expresa por número de ASVs y lecturas asociadas en cada
rango/umbral, con cobertura adicional por muestra. Las ASVs con nombre
`uncultured` cuentan como etiquetas asignadas, pero no como géneros resueltos.
Las figuras de cobertura se exportan en PDF, SVG y PNG de 300 dpi.

## Reproducción y validación

```bash
python3 scripts/environment/run_in_environment.py snakemake --snakefile workflow/Snakefile --cores 2 pilot_taxonomy
python3 scripts/environment/run_in_environment.py python scripts/taxonomy/validate_taxonomy.py
python3 scripts/environment/run_in_environment.py Rscript tests/taxonomy/test_helpers.R
```

La regla de descarga reconstruye y verifica la referencia a partir del DOI,
URL, bytes y MD5 fijados. La base se protege y se separa de la regeneración
del registro de validación: cambiar parámetros o código no reemplaza la base.
El ejecutor verifica previamente los resultados
DADA2 y compara hashes de todas sus entradas antes y después de clasificar.
El registro conserva commit, estado inicial de Git, semilla, parámetros,
versiones de R/paquetes y SHA-256 de entradas y artefactos. `SUCCESS` se crea
solo al terminar; una nueva ejecución invalida el marcador anterior.

El validador Python independiente vuelve a calcular máscaras de bootstrap,
marcas, pertenencia de ASVs y cobertura global/por muestra a partir de las
tablas R y conteos DADA2. También comprueba los hashes de entradas y salidas.
Las pruebas con fixtures sintéticos se limitan a parsing y conservación;
no aportan evidencia científica.
