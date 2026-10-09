# Métodos de separación bacteriana y diversidad del piloto

## Alcance

Esta etapa continúa el piloto de seis runs de PRJNA492720, con cuatro lotes
identificados. Las muestras y sus estratos permanecen separados; no se
construye una trayectoria completa de cada lote ni se estiman efectos
temporales. Los parámetros se registraron en `f9132b5` antes de calcular
diversidad, en `config/diversity.yaml` y `protocol/analysis_decisions.md`.

## Separación auditable

Se consume la llamada taxonómica primaria con bootstrap 80, previamente
validada. Se separan etiquetas exactas Chloroplast/Mitochondria en cualquier
rango y Kingdom conocido distinto de Bacteria. La coincidencia por rango
ignora mayúsculas sólo para orgánulos; no busca subcadenas. Kingdom desconocido
se conserva y marca, al igual que etiquetas inferiores sin asignación. No
hay filtro por abundancia, prevalencia, etapa o significación. Una muestra
vacía causa fallo explícito, no una exclusión silenciosa.

El [tutorial oficial de QIIME 2](https://docs.qiime2.org/2024.10/tutorials/filtering/)
documenta la separación taxonómica de orgánulos. Aquí se implementa sobre
rangos SILVA originales sin instalar QIIME 2. No equivale a identificar
contaminación o especie huésped; no hay controles negativos verificados.

Las tablas originales DADA2/taxonómicas se preservan. Los conteos derivados
quedan en `data/processed/pilot/bacterial/`, conservando IDs, secuencias,
etiquetas y metadata. `results/filtering/pilot/asv_filter_log.tsv` contiene
todas las ASVs, incluidas las retenidas, con razón y lecturas originales.
`sample_retention.tsv` comprueba para cada muestra:

```text
lecturas originales = lecturas retenidas + lecturas separadas
```

El validador independiente vuelve a decidir a partir de la taxonomía original
y comprueba cada celda retenida, pertenencia, secuencia, etiqueta y campo de
metadata. No usa el resultado del productor para decidir qué se debió excluir.

## Alfa diversidad descriptiva

Se usan los conteos retenidos sin rarefacción. Para una muestra, `p_i` es la
proporción de su ASV i. Se exportan riqueza observada (ASVs con conteo >0),
Shannon `H = −Σ p_i ln(p_i)`, Gini-Simpson `1−Σ p_i²`, inverso de Simpson
`1/Σ p_i²` y diversidad efectiva `exp(H)`, según las definiciones de
[vegan](https://vegandevs.github.io/vegan/reference/diversity.html).

Son índices de bibliotecas observadas. Profundidad desigual y detección de
ASVs raras afectan su comparación; no estiman riqueza total ni aportan una
prueba de diferencias entre etapas. Los gráficos incluyen profundidad y
muestras individuales. Se conserva explícitamente el índice Simpson utilizado
para evitar confundir concentración, complemento e inverso.

## Beta diversidad y PCA

Para CLR se suma un pseudoconteo de 1 a **todos** los conteos; sensibilidad 0,5
sobre las mismas ASVs. Se calcula `z_i = ln(x_i+c) − media_j ln(x_j+c)` y la
distancia euclídea entre vectores CLR (Aitchison de las composiciones
sustituidas), conforme a la [definición de CLR](https://scikit.bio/docs/dev/generated/skbio.stats.composition.clr.html).
No se añaden lecturas a los conteos originales ni a los índices alfa.
La sustitución depende de escala/profundidad y no identifica ceros
estructurales; se reportan ambos valores sin elegir por separación visual.

El diagnóstico secundario Bray-Curtis usa proporciones sin pseudoconteo.
PCA centra columnas CLR y no estandariza varianzas. Se exportan hasta
`min(n−1, número de ASVs)` componentes, varianzas y fracciones explicadas.
El signo se fija mediante la carga de mayor valor absoluto; su orientación
no implica una dirección biológica. La figura muestra las primeras dos
componentes y ambos pseudoconteos, sin unir muestras en una trayectoria.

Python comprueba fórmulas alfa/CLR/distancias, identificadores, finitud,
centrado y ortogonalidad de scores PCA, orden de varianzas y conservación de
las distancias al usar todas las componentes. Se registran configuración,
semilla 20260819, versiones instaladas y hashes antes/después. Los métodos
actuales son deterministas; no se simula ni se submuestrea.

## Ejecución

```bash
python3 scripts/environment/run_in_environment.py snakemake --snakefile workflow/Snakefile --cores 2 pilot_diversity
python3 scripts/environment/run_in_environment.py python scripts/filtering/validate_bacterial_table.py
python3 scripts/environment/run_in_environment.py python scripts/diversity/validate_diversity.py
```

El target `pilot_bacterial_table` detiene el workflow en la separación validada.
El target por defecto `all` alcanza diversidad y sus validaciones. Agregar
estas reglas cambia el archivo principal de workflow, registrado como entrada
taxonómica; por eso se reejecuta esa clasificación con los mismos parámetros.
Se comparan sus siete tablas científicas contra el hito anterior por SHA-256.

No se ejecutan PERMANOVA, abundancia diferencial, core ni meta-análisis en el
piloto técnico. El procesamiento completo por estudio, tratamiento de
submuestras/profundidad y diagnóstico de dispersión siguen siendo requisitos
antes de interpretar sucesión. Los fixtures sintéticos sólo prueban contratos
y fórmulas; no aparecen como resultados del estudio.

## Resultados observados y validación

La ejecución del 2026-10-09 sobre `5221b0f` conserva seis muestras y separa
cinco ASVs etiquetadas como orgánulos, con 32.868 lecturas. La tabla derivada
contiene 88 ASVs y 123.693 lecturas, sin exclusiones de muestras ni Kingdom
desconocido. La separación valida 21 entradas y nueve artefactos. Las siete
tablas taxonómicas reproducen por SHA-256 las del hito anterior.

La diversidad valida 20 entradas y 17 artefactos. Profundidades: 9.250–36.508
lecturas; riqueza observada: 2–46 ASVs; diversidad efectiva de Shannon:
1,013–6,072. PC1+PC2 explican 68,92 % con pseudoconteo 1 y 68,87 % con 0,5.
Se exportan las cinco componentes para verificar las distancias completas;
estas fracciones no miden evidencia de sucesión ni significación.

Las versiones son R 4.4.3, vegan 2.6.8 y yaml 2.3.10, sin warnings R. Los
ajustes de presentación reservan espacio para etiquetas y títulos y usan
la misma escala geométrica por unidad en los dos ejes PCA. Las seis tablas
científicas y la tabla de versiones reproducen sus hashes en la ejecución
final sobre `342dd4e`, registrada en `results/diversity/pilot/provenance.json`.
Las figuras finales se revisaron visualmente y el dry-run integrado quedó
sin trabajos pendientes.
La suite completa aprobó 125 pruebas y las 12 pruebas específicas de esta
etapa se repitieron tras los ajustes de figuras.

El estado Git inicial se conserva íntegro en cada procedencia. Snakemake
retira salidas previas antes de regenerarlas; las etapas posteriores observan
también resultados recién producidos y aún sin commit. `git_dirty=true` no se
oculta: los commits de código y los hashes de entradas/salidas permiten
reconstruir qué se ejecutó. No había cambios sin commit de código o
configuración al lanzar las ejecuciones científicas.
