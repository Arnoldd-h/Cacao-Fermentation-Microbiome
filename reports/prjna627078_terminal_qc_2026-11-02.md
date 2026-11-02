# QC terminal de PRJNA627078 — 2026-11-02

La segunda pasada conserva las 60 bibliotecas y los 5.689.175 pares del QC
basal. Retira coincidencias completas del primer opuesto en el extremo 3',
sin filtrar lecturas ni modificar las calidades de las bases conservadas.
No produce todavía ASVs ni evidencia de sucesión compartida entre estudios.
Las cuatro cajas siguen siendo las unidades independientes; las 60 bibliotecas
no representan 60 fermentaciones independientes.

## Método y trazabilidad

Código y decisión registrados antes de la ejecución: `c77233d6a9acbd2ad45f15cf4433db30cdf96a82`.
Configuración: [`prjna627078_terminal_qc.yaml`](../config/prjna627078_terminal_qc.yaml).
El baseline 5' de `full` se conserva; salidas nuevas en `terminal_3prime`.

Cutadapt 5.0 aplica una sola pasada, error 0,1, sin indels: complemento inverso
del primer reverse en R1 (21 nt) y del forward en R2 (17 nt), completos y
anclados con `$`. Se conservan coincidencias internas y parciales. El anclaje
exige el sufijo completo y no admite `min_overlap` explícito en esta versión.
[Documentación de Cutadapt 5.0](https://cutadapt.readthedocs.io/en/v5.0/guide.html#anchored-3-adapters).

La auditoría compara los cuatro streams baseline/post en orden: identidad,
cantidad de pares, prefijos de secuencia/calidad y recorte de exactamente
0 o 21/17 nt. Contrasta reads y bases con Cutadapt; después criba por separado
el primer opuesto completo terminal y en cualquier posición. Cada run conserva
hashes del método, entradas, salidas y cinco informes de ejecución.

FastQC 0.12.1 se repite sobre los 120 archivos. Sus intervalos se acompañan de
cobertura exacta calculada desde las longitudes: reads que alcanzan inicio/final
y observaciones de bases. Esto no reconstruye ni pondera de nuevo sus cuantiles.
Los `NaN` nativos permanecen como estadísticas no disponibles, sin imputación.

## Resultados observados

| Indicador | R1 | R2 |
|---|---:|---:|
| Lecturas examinadas por mate | 5.689.175 | 5.689.175 |
| Lecturas con sufijo retirado | 92.518 | 76.111 |
| Bases retiradas | 1.942.878 | 1.293.887 |
| Coincidencias terminales posteriores | 0 | 0 |
| Coincidencias completas posteriores en cualquier posición | 8.961 | 7.513 |
| Lecturas posteriores de 148 nt | 85.356 | 72.914 |
| Rango posterior de longitud, nt | 3–251 | 3–251 |

Los recuentos por mate no se suman como pares únicos. Los ceros terminales
corresponden al primer opuesto y a la tolerancia registrada; no significan
ausencia de cualquier primer, coincidencia parcial o secuencia técnica.
El cribado anclado y el diagnóstico ilustrativo previo seleccionan alineamientos
distintos y no tienen por qué producir idénticos recuentos terminales.

El mínimo candidato de DADA2 sigue en 150 nt, sin aplicarse aquí:

| Pares con algún mate menor de 150 nt | Recuento |
|---|---:|
| Baseline 5' | 6.088 |
| Después del recorte terminal | 95.898 |
| Nuevos cruces del umbral | 89.810 |

Quedan 5.593.277 pares con ambos mates de al menos 150 nt (98,3144 % de entrada).
Esto no estima la retención posterior por calidad, ensamblaje o quimeras.
SRR11579032 conserva sus 877 pares; 873 cumplen ese criterio de longitud.
La menor proporción es SRR11579065: 178.642/192.977 (92,5717 %).
Ninguna biblioteca queda sin pares que cumplan el criterio. No se baja el umbral
para recuperar el pico de 148 nt ni se asigna identidad biológica a ese pico.

## Calidad y cobertura

Los flags por calidad de base permanecen iguales al baseline 5':

| Mate | Pass | Warn | Fail |
|---|---:|---:|---:|
| R1 | 58 | 1 | 1 |
| R2 | 0 | 0 | 60 |

Las tablas contienen 13.920 intervalos de calidad (120 archivos por dos etapas).
Todos los últimos intervalos abarcan posiciones 250–251. Después del recorte,
la fracción de reads que alcanza la posición 250 tiene mediana entre bibliotecas
de 1,8916 % en R1 y 3,0500 % en R2; en R2 varía entre 2,7171 y 7,6397 %.
Esas son estadísticas de cobertura entre archivos, no cuantiles de calidad
combinados entre muestras.

En SRR11579065/R2, el último intervalo conserva mediana Q10 y recibe 7.949
lecturas (4,1191 %), frente a 7.951 en el baseline. Los percentiles ausentes
del archivo de menor profundidad no se interpretan como Q0. En toda la cohorte,
181.291 R2 posteriores superan 230 nt. Las colas corresponden a subconjuntos
de longitudes distintas; sus flags no justifican excluir bibliotecas ni
describen por sí solos la calidad de todos los pares.

## Validación y relación con el objetivo

La suite completa aprueba 217 pruebas, incluidas 24 nuevas de método terminal,
conservación, corrupción, sincronía, coincidencias internas/parciales,
degeneración, tolerancia y cobertura incompleta. Un cálculo independiente
reconcilia los 60 runs, 300 informes y los agregados anteriores.
La procedencia registra 861 entradas y seis artefactos con SHA-256, código
productor `c77233d` y checkout limpio al iniciar el cierre.

El avance reduce una fuente técnica de diferencias entre estudios y cuantifica
su efecto sobre la elegibilidad de lecturas antes de interpretar microbiología.
Persisten coincidencias no terminales y alertas de calidad. DADA2 mantiene
`processing_ceiling=qc`: se debe registrar la ejecución que consuma
`terminal_3prime`, revisar los filtros candidatos (minLen 150, maxEE 2/2,
truncQ 2, sin truncamiento fijo) y verificar convergencia, ensamblaje y quimeras.
La ruta independiente debe adaptar las tablas QC al runner existente y validar
la presencia de las 60 bibliotecas; no cambiar el scope basal ni reutilizar
el modelo colombiano de estratos espaciales desconocidos en México.

```bash
python3 scripts/environment/run_in_environment.py snakemake --snakefile workflow/terminal_qc.Snakefile --cores 6 --resources mem_mb=3500 --keep-going
python3 scripts/environment/run_in_environment.py python scripts/terminal_qc/run.py validate --config config/prjna627078_terminal_qc.yaml
```

Resultados: [`summary.json`](../results/qc/mota_gutierrez_2021_mexico/terminal_3prime/summary.json),
[`run_summary.tsv`](../results/qc/mota_gutierrez_2021_mexico/terminal_3prime/run_summary.tsv),
[`read_length_counts.tsv`](../results/qc/mota_gutierrez_2021_mexico/terminal_3prime/read_length_counts.tsv),
[`quality_coverage.tsv`](../results/qc/mota_gutierrez_2021_mexico/terminal_3prime/quality_coverage.tsv).
