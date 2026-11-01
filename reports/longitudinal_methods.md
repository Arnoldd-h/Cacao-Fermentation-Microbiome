# Modelo longitudinal dentro del primer estudio

Plan registrado el 2026-10-31, después de conocer QC y resultados descriptivos
y **antes del primer ajuste longitudinal**. Es un análisis exploratorio de un
estudio; no una prerregistración prospectiva externa ni un meta-análisis.
No cambia las hipótesis del protocolo ni selecciona Acetobacter por el patrón
ya observado. `config/longitudinal.yaml` contiene todas las decisiones.

## Unidad y población analizada

Se utilizan las ASVs bacterianas retenidas del estudio completo, sin unirlas
con el piloto ni con otro estudio. Las cuatro fermentaciones son las unidades
independientes. Hay 93 muestras, 47 combinaciones lote/hora y dos estratos
esperados por observación. BactAnII108z1 sólo dispone del estrato 1: la ruta
pareada no pudo procesar el run del estrato 2. No se inventa una pareja.

El principal exige los dos estratos y utiliza 46 observaciones/92 muestras.
La muestra restante permanece en los datos fuente y en la sensibilidad que
usa todos los estratos disponibles. También se ajusta cada estrato por separado,
siempre con cuatro lotes; sus dos resultados no constituyen ocho réplicas.
No se excluyen otros puntos por abundancia, profundidad, QC o resultado.
Cada pertenencia/exclusión por escenario tendrá una fila con motivo auditable.

## Transformación y ajuste

Para cada muestra y cada una de las D ASVs se calcula
`z_j = log(count_j + c) - mean_k(log(count_k + c))`. El principal usa c=1;
c=0,5 es sensibilidad. Son los pseudoconteos ya utilizados en diversidad.
Se conservan todas las ASVs, incluidos ceros originales, sin rarefacción ni
selección por prevalencia. La geometría CLR/Aitchison refleja relaciones entre
componentes, no crecimiento absoluto; su fundamento composicional se describe
en [Gloor et al.](https://doi.org/10.3389/fmicb.2017.02224).

Dentro de cada lote/hora se promedian los vectores CLR con igual peso por
estrato. Es el centro de Aitchison de composiciones suavizadas, no una suma
de lecturas ni una ponderación por profundidad. El pseudoconteo sigue siendo
sensible a la profundidad; la sensibilidad de ceros y los estratos lo hacen
visible, no lo eliminan. Las riquezas alfa no se agregan ni modelan en esta fase.

Se ajusta, por separado para cada lote, el modelo multivariado lineal con
intercepto `z_bt = a_b + x_bt * beta_b + e_bt`, con igual peso por hora.
Cada coordenada comparte el mismo diseño OLS. La estimación combinada es
`mean_b(beta_b)`, con igual peso por lote, independientemente de sus números
de tiempos. Se exportan todos los coeficientes ASV; no se les asignan valores p
individuales ni se usan como prueba de abundancia diferencial por taxón.

El principal usa el tiempo relativo de la metadata, cuya duración procede del
último tiempo observado (108–144 h). La sensibilidad absoluta usa horas/24.
No se interpreta el último tiempo como un estado biológico equivalente y las
magnitudes por duración observada y por 24 horas tienen unidades diferentes.
El modelo lineal representa una dirección media; no prueba monotonicidad ni
descarta cambios no lineales o transitorios. Se reportan ajuste y residuos por
lote para mostrar sus límites, sin seleccionar otro modelo por su resultado.

## Contraste global e incertidumbre

La estadística es `T = norm(mean_b(beta_b))` en geometría CLR. Se enumeran
todos los cambios de signo de los vectores completos de pendiente por lote,
sin permutar tiempos ni estratos dentro de las series. El valor p es la
fracción inclusiva de estadísticas al menos tan grandes como T; no hay
muestreo Monte Carlo ni corrección +1 en una enumeración exhaustiva.

La validez condicional exige lotes independientes y simetría central de los
vectores de pendiente alrededor de cero bajo la nula. No se supone
intercambiabilidad entre horas autocorrelacionadas. Estos supuestos no pueden
comprobarse con cuatro lotes y no convierten el estudio observacional en un
experimento aleatorizado. El principio de signos independientes/simétricos
se apoya en [Winkler et al.](https://doi.org/10.1016/j.neuroimage.2014.01.060);
la estadística vectorial concreta es la decisión explícita de este proyecto.

Con cuatro lotes hay 2^4=16 patrones. Como T es invariante al signo global,
el mínimo alcanzable es 2/16=0,125. No es posible rechazar a alfa 0,05 con
este diseño y contraste. Un valor p mayor de 0,05 no demostrará ausencia de
sucesión; la resolución disponible no permite esa conclusión. No se ampliará
artificialmente el número de réplicas para obtener significación.

Se reportan magnitudes por lote, cosenos entre las seis parejas de direcciones,
dispersión cuadrática de pendientes respecto a su media y cuatro resultados
leave-one-batch-out. Se conservan vectores de efecto y su variación, sin IC
paramétricos de alta precisión ni Q/I² de un meta-análisis inexistente.
El contraste de signos admite dependencia temporal interna bajo la hipótesis
de simetría del vector por lote; esa hipótesis sigue siendo una limitación.

Los seis escenarios globales forman una familia. Se publica p bruto y FDR
Benjamini–Yekutieli (BY), conservadora ante dependencia entre sensibilidades,
según la [documentación primaria de R](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/p.adjust.html).
No se elige el escenario con menor p. Las figuras de distancia desde 0 h son
descriptivas: alejarse de un origen no demuestra compartir dirección.

## Escenarios fijados

| Escenario | Tiempo | c | Estratos |
|---|---|---:|---|
| Principal | Relativo | 1 | Ambos completos |
| Ceros | Relativo | 0,5 | Ambos completos |
| Horas absolutas | Por 24 h | 1 | Ambos completos |
| Disponibles | Relativo | 1 | Todos los disponibles |
| Estrato 1 | Relativo | 1 | Sólo 1 |
| Estrato 2 | Relativo | 1 | Sólo 2 |

Se exige como mínimo 0 h y tres horas distintas por lote/escenario. Metadata
faltante, duplicación lote/hora/estrato, mezcla de estudios, un lote perdido o
un diseño sin variación temporal producen fallo; no se imputan valores ni
se relajan requisitos después de observar efectos. El piloto técnico queda
fuera del análisis inferencial. Semilla 20260819 registrada, aunque las
operaciones y la enumeración de esta fase son deterministas.

## Verificación y alcance

Las entradas deben pasar la validación de separación bacteriana y sus hashes.
Antes de analizar los datos reales se usarán fixtures sintéticos con solución
analítica para OLS, centrado CLR, promedio de estratos, pesos entre lotes,
signos exhaustivos y FDR. Una comprobación independiente en R con `lm.fit` y
`p.adjust` contrastará las fórmulas y la pertenencia de observaciones.
Salidas y figuras conservarán configuración, versión Git y SHA-256.

Esta fase estima estructura temporal dentro de cuatro lotes del primer estudio.
No identifica efectos independientes de país, temporada, región 16S o tecnología;
no define un core ni prueba sucesión conservada entre estudios. Esas hipótesis
requieren procesar los otros estudios por separado y armonizar taxonomía.
