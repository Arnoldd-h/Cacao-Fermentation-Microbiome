# Figuras temporales descriptivas por lote y estrato

La política se registró en `bfd48e5`, en `config/temporal_reporting.yaml` y
`protocol/analysis_decisions.md`. Se representan tres columnas de diversidad
alfa previamente validada: ASVs observadas, diversidad efectiva de Shannon y
lecturas retenidas. No se recalculan métricas ni se modifican valores, ceros,
etapas o tiempos. Lecturas retenidas indica profundidad de biblioteca, no
abundancia celular absoluta.

Cada lote tiene paneles separados y cada estrato conserva su identidad y color.
El eje X utiliza las horas originales y las escalas Y son comunes entre lotes
para cada métrica. No se agregan submuestras, no hay rarefacción ni pruebas
estadísticas. La profundidad acompaña a la riqueza para hacer visible su
limitación de comparación.

En el perfil completo se unen puntos únicamente dentro del mismo lote/estrato
y cuando son consecutivos entre los tiempos observados del lote. Si falta una
hora en un estrato, la línea se interrumpe; no se dibuja un punto imputado ni
se conecta el salto. Las líneas son guías entre observaciones, no modelos.
El piloto se muestra con puntos aislados porque no reconstruye series completas.
Duplicados de muestra, run o lote/hora/estrato provocan fallo explícito.

La tabla `plot_data.tsv` conserva exactamente los identificadores, metadata y
valores utilizados, ordenados por lote/hora/estrato. El productor valida de
nuevo la etapa fuente y registra hashes de entradas, código, configuración,
salidas y sesión R. Las figuras se exportan en PDF/SVG/PNG de 300 dpi; el
validador comprueba integridad, pertenencia exacta de puntos y dimensiones PNG.
Se requiere revisión visual adicional antes de considerar las figuras finales.

Los tests sintéticos comprueban conservación de valores, separación de estratos,
identidades duplicadas, tiempos no finitos, ausencia de agregación y saltos en
las líneas R. No se presentan como evidencia biológica.

Después de completar la diversidad del perfil elegido:

```bash
python3 scripts/environment/run_in_environment.py python scripts/reporting/plot_temporal_diversity.py --analysis full_study
python3 scripts/environment/run_in_environment.py python scripts/reporting/plot_temporal_diversity.py --analysis full_study --validate-only
```

`--analysis pilot` utiliza la misma política, con líneas desactivadas. No se
añaden dependencias: se usan Python estándar, R base y `yaml` ya declarados.
Estos gráficos no demuestran sucesión conservada entre estudios ni definen un core.

Primera ejecución piloto: seis muestras, cuatro lotes y seis series parciales,
sobre `afb7705`. Se validaron 16 entradas y ocho artefactos; los valores de
`plot_data.tsv` coinciden exactamente con la diversidad fuente. La revisión
visual del PNG comprobó títulos, ejes, leyenda, colores y puntos sin conexiones.
Se corrigieron los márgenes antes de aceptar esta exportación.

El perfil completo contiene 93 muestras, cuatro lotes y ocho series de
lote/estrato. Se validan 16 entradas y ocho artefactos. La revisión visual
del PNG comprueba títulos, ejes, leyenda, separación de lotes y correspondencia
de las trayectorias con sus horas observadas. Se mantienen huecos sin imputar
y la profundidad junto a la riqueza. La fuente exacta queda en la procedencia
de la ejecución vigente; las figuras son descriptivas.
