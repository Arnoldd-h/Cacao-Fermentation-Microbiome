# Primer modelo longitudinal: PRJNA492720

Ejecución real el 2026-10-31 local, fuente
`7b6d5a6570034f3ba523256254860cca32afde2a`, checkout limpio al iniciar.
El plan se fijó previamente en `c97345c`, después de la descripción inicial.
La salida actual procede de la renovación operativa `6415bfd`, sin cambiar
el análisis ni las tablas científicas de la primera ejecución.
Métodos y supuestos: [longitudinal_methods.md](longitudinal_methods.md).

## Diseño y validación

El principal analiza 764 ASVs, 46 observaciones lote/hora y 92 muestras de
cuatro lotes independientes. La única muestra omitida del principal es
BactAnII108z1, porque falta su otro estrato en la ruta pareada; permanece en
la sensibilidad de estratos disponibles (47 observaciones/93 muestras).
Los escenarios por estrato usan 47 y 46 muestras, con los mismos cuatro lotes.
La tabla de pertenencia registra las 558 decisiones muestra/escenario.

El primer recorrido cerró tres trabajos. Después se unificó la propiedad del
archivo de validación en el productor, evitando un aviso del scheduler y una
comprobación duplicada. La renovación sobre `6415bfd` cerró dos trabajos,
sin regenerar QC, DADA2, taxonomía, separación ni diversidad. Las diez tablas
científicas, resumen, comprobación R y figuras PNG/SVG son idénticos por SHA-256
a la primera ejecución completa archivada. La única entrada modificada fue
`workflow/longitudinal.Snakefile`; no cambió un parámetro ni se eligió una
corrida por resultado. La nueva procedencia conserva `git_dirty=true` por las
notas/resultados pendientes de commit, con todas las entradas analíticas
fijadas por hashes. La validación comprueba 27 entradas y
19 artefactos. R reconstruyó los centros desde los conteos y contrastó 36.672
coeficientes de pendiente/intercepto con QR, además de membresía, residuos,
estadísticas de signos y FDR BY en seis escenarios. No hubo warnings R en
la ejecución. La figura PDF/SVG/PNG de 300 dpi fue revisada visualmente.

## Resultados observados

La norma de la pendiente media es 12,152 unidades CLR por duración observada.
La media de normas individuales es 16,679 y la dispersión RMS entre pendientes
es 11,976: las trayectorias presentan variación entre lotes. Estas magnitudes
dependen del conjunto fijo de 764 componentes y no se comparan directamente
con ASVs de otro estudio.

| Lote | Tiempos analizados | Norma de pendiente | R² multivariado lineal |
|---|---:|---:|---:|
| An::I | 12 | 14,135 | 0,1205 |
| An::II | 12 | 13,710 | 0,1182 |
| Sa::I | 10 | 16,176 | 0,1644 |
| Sa::II | 12 | 22,693 | 0,1284 |

Los seis cosenos entre direcciones de pendiente son positivos, entre 0,249
y 0,584. Esto describe alineación parcial de la dirección temporal dentro
del estudio. El ajuste lineal recoge sólo el 11,82–16,44 % de la variación
CLR por lote; no resume todos los cambios ni prueba monotonicidad. Las
distancias observadas desde 0 h tienen transitorios y no equivalen por sí
solas a una dirección común.

Al omitir un lote a la vez, la norma de la pendiente media varía entre 11,992
y 13,226 y su coseno con la dirección completa entre 0,898 y 0,973. Es una
sensibilidad interna con tres lotes restantes, no cuatro estudios nuevos
ni un intervalo de confianza poblacional.

| Escenario | Observaciones / muestras | Norma de pendiente media | p exacto | q BY |
|---|---:|---:|---:|---:|
| Principal relativo, c=1 | 46 / 92 | 12,152 | 0,125 | 0,30625 |
| Relativo, c=0,5 | 46 / 92 | 13,386 | 0,125 | 0,30625 |
| Horas/24, c=1 | 46 / 92 | 2,294 por 24 h | 0,125 | 0,30625 |
| Estratos disponibles | 47 / 93 | 12,012 | 0,125 | 0,30625 |
| Estrato 1 | 47 / 47 | 13,805 | 0,125 | 0,30625 |
| Estrato 2 | 46 / 46 | 13,459 | 0,125 | 0,30625 |

Los seis valores p alcanzan el mínimo posible de este contraste con cuatro
lotes: 2/16=0,125. Ninguno rechaza a alfa 0,05. La resolución impide alcanzar
ese umbral incluso con direcciones perfectamente alineadas; por ello la no
significación no demuestra ausencia de cambio temporal. Las condiciones de
independencia y simetría central bajo la nula tampoco se verifican con cuatro
unidades. No se cambió el contraste, el conjunto de ASVs ni las sensibilidades
para buscar un valor menor.

## Relación con el objetivo original

El proyecto ya dispone de una estimación explícita y auditable de dirección
temporal dentro del primer estudio, más allá de la comparación descriptiva
de abundancias o de una figura PCA. La alineación parcial entre cuatro lotes
es compatible con la hipótesis de un componente de sucesión compartido,
pero no la demuestra entre estudios o geografías. La heterogeneidad y el
modesto ajuste lineal permanecen visibles.

La observación anterior de Acetobacter no se utilizó para seleccionar el
contraste y no recibe aquí una prueba individual de abundancia diferencial.
El core temporal, los efectos por taxón/etapa y la conservación cross-study
siguen pendientes de datos independientes y modelos registrados apropiados.
No se reemplazará este resultado por un modelo seleccionado por significación.

Las bibliotecas originales varían de 986 a 69.751 lecturas retenidas. El peso
igual por estrato evita que la biblioteca mayor domine directamente el centro,
pero no elimina la dependencia del pseudoconteo respecto a la profundidad.
La sensibilidad c=0,5 no constituye una corrección completa de ese problema.
Las pendientes no se interpretan como crecimiento absoluto ni efectos causales.

## Salidas reconstruibles

Todo se genera con `workflow/longitudinal.Snakefile`, después del workflow
del estudio completo. Los archivos están en
`results/longitudinal/pacheco_montealegre_2020_colombia/full/`:

- `global_tests.tsv`, `batch_models.tsv`, `slopes.tsv` e `intercepts.tsv`.
- `sample_membership.tsv` y `observations.tsv`: pertenencia y estratos sin imputar.
- `direction_cosines.tsv`, `leave_one_batch_out.tsv` y `sign_statistics.tsv`.
- `trajectory_plot_data.tsv` y `community_trajectories.pdf/.svg/.png`.
- `independent_validation.tsv`, versiones, configuración y procedencias SHA-256.

Validación del desarrollo: suite de 181 pruebas (152,629 s), más dos pruebas
del nuevo workflow (10,738 s en su repetición final), 183 pruebas distintas
aprobadas. Las dos últimas
verifican aislamiento de etapas previas e invalidación por contenido de
configuración incluso cuando su timestamp es anterior a las salidas.
La comparación de renovación está en
`results/reproducibility/pacheco_montealegre_2020_colombia/longitudinal_workflow_renewal.json`.
El resumen verificable de pruebas/validadores/workflows está en
`results/reproducibility/pacheco_montealegre_2020_colombia/longitudinal_validation_summary.json`.
La entrega corresponde a [PR #8](https://github.com/Arnoldd-h/Cacao-Fermentation-Microbiome/pull/8)
y al tag `v0.5.3-longitudinal`.
