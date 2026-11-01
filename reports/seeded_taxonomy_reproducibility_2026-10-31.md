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

La evaluación real debe comparar cronológicamente las siete tablas científicas
en procesos independientes, contar las celdas afectadas frente al histórico y
renovar derivados con la corrida definida, sin seleccionar por cobertura o
significación. Su ejecución y resultados se añadirán a este registro.
