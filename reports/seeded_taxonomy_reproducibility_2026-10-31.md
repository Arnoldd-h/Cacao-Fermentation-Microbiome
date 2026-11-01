# Control explícito del desempate taxonómico

Fecha local: 2026-10-31. La decisión se registra antes de reclasificar datos
reales en `protocol/analysis_decisions.md` y `config/taxonomy.yaml`.

El código original de [DADA2 1.34.0](https://github.com/bioc/dada2/blob/5d453a416608f75f0bf50bcc74bbf0c1851e7f9f/src/taxonomy.cpp)
resuelve empates con `std::random_device`. El adaptador explícito del proyecto
conserva el clasificador y cambia esa fuente por un generador privado por
consulta, inicializado desde R después del vector bootstrap. Se preservan
SILVA 138.2, la semilla 20260819, 100 réplicas y umbrales 80/50. El paquete
instalado no cambia; fuente LGPL-2 original y modificada quedan versionadas.

## Validación previa

Tres procesos independientes clasifican las mismas 300 consultas sintéticas
contra dos linajes con puntuaciones deliberadamente empatadas, con uno, dos
y dos hilos. Las tablas y soportes coinciden exactamente. Ambos géneros son
seleccionados y el soporte refleja la ambigüedad; no se fija la primera etiqueta.
Un caso sin linajes empatados concuerda con la entrada nativa instalada.
Estos fixtures no se presentan como evidencia biológica.

Las tablas completas originales del piloto y estudio completo están archivadas
en `results/tmp/taxonomy_native_baseline`; el registro versionado
`results/reproducibility/pacheco_montealegre_2020_colombia/taxonomy_native_baseline.json`
conserva hashes de todos sus archivos y el commit original `cb9eca0`, desde el
que también se recuperan las tablas históricas. No se sobrescribe ese archivo.

## Alcance de la reproducción

Se exige el mismo orden de consultas y referencias, semilla, configuración y
entorno registrado. La prueba entre hilos evalúa independencia del reparto
de trabajo; no promete igualdad entre compiladores/plataformas ni al reordenar
secuencias. Cada proceso recompila explícitamente para evitar binarios de
caché ajenos al compilador registrado. Las versiones y flags quedan en la tabla
de software y las fuentes en los hashes de entrada.

## Ejecución real y resultados

Fuente científica: `e008f2b59c17b46915e5dcca6b6a4da0231befb8`. Se ejecutó primero
el workflow completo con dos hilos taxonómicos y después un proceso independiente
con un hilo; se repitió esa secuencia en el piloto. Cada workflow cerró siete
trabajos sin repetir QC ni DADA2. El workflow renueva automáticamente separación
y diversidad antes de la auditoría final; no se interpretaron ni publicaron
esos derivados antes de comprobar la comparación. Esta es una precisión del
orden operativo previsto; no altera la decisión analítica ni selecciona corridas.

Las siete tablas científicas de cada alcance son idénticas por SHA-256 entre
dos y un hilo, con iguales entradas/configuración comprobadas. Las ejecuciones
independientes pasan el validador de 26 entradas y 15 artefactos, al igual que
las tablas actuales. Las fechas UTC e hilos se conservan en los dos registros
`taxonomy_seeded_*_repeat.json`. El archivo completo de cada repetición queda
en `results/tmp/taxonomy_seeded/`; las tablas repetidas se recuperan también
de las salidas actuales, cuyos hashes coinciden exactamente.

Frente al histórico completo, sólo cambian dos celdas bootstrap en una ASV:
`pacheco_montealegre_2020_colombia__ASV000100`, cinco lecturas en BactSaI012z1.
Familia y género pasan de soporte 48 a 51. A 50 aparecen Yersiniaceae/Samsonia;
ambas etiquetas siguen vacías en la tabla primaria a 80. La clasificación sin
máscara, la primaria y el screening son idénticos. Cambian además dos celdas
de sensibilidad, cuatro de cobertura global y dos de cobertura por muestra,
todas trazables a esa ASV. No se interpreta la sensibilidad como identificación
primaria ni se atribuye relevancia biológica a cinco lecturas.

| Alcance | Género a 80 | Género a 50 | Comparación con histórico |
|---|---|---|---|
| Completo, 782 ASVs | 475 / 1.855.013 lecturas | 649 / 2.051.014 lecturas | Primaria idéntica; sensibilidad cambia una ASV/cinco lecturas |
| Piloto, 93 ASVs | 64 / 99.481 lecturas | 79 / 104.974 lecturas | Siete tablas científicas idénticas |

El registro `taxonomy_seeded_derived_comparison.json` conserva 44 TSV idénticos
de DADA2, tablas bacterianas y diversidad en ambos alcances. No cambian las
764 ASVs/2.096.980 lecturas retenidas del estudio completo ni sus 93 muestras.
`taxonomy_seeded_reporting_comparison.json` conserva ocho TSV/PNG de composición
y figuras temporales. Los PDF renovados pueden diferir por metadatos del
dispositivo; los PNG científicos coinciden. No hay warnings R de clasificación
o diversidad. Las procedencias de las etapas renovadas describen su nueva
clasificación de origen; no se editan resultados manualmente.

Se conserva la primera corrida controlada como salida del workflow y se
verifica con la siguiente; no se selecciona una réplica por cobertura o
significación. La coincidencia con el total de una antigua corrida nativa
no fue criterio de aceptación. Esto cierra el problema observado dentro del
contrato de entorno/orden, sin garantizar reproducción universal de DADA2.

## Cierre de validación

La suite completa pasa 166 pruebas en Linux (140,517 s). Se añade después una
regresión de integridad de auditoría: tablas iguales con entradas distintas
deben fallar; pasa junto a los otros dos casos de comparación. Las pruebas
del nuevo clasificador incluyen cambio de semilla, conservación del namespace
instalado, tres procesos y uno/dos hilos. La referencia fue revalidada sin
descarga: 452.055 secuencias, mismos bytes/MD5. La nueva figura de cobertura
se revisó visualmente; los PNG de diversidad/tiempo conservan los hashes
previamente revisados. Los dos workflows quedan sin trabajos pendientes.
