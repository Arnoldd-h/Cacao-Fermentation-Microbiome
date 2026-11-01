# Adaptador de desempate reproducible

`taxonomy_upstream.cpp` es una copia sin modificar de DADA2 1.34.0,
Bioconductor 3.20, commit `5d453a416608f75f0bf50bcc74bbf0c1851e7f9f`:
https://github.com/bioc/dada2/blob/5d453a416608f75f0bf50bcc74bbf0c1851e7f9f/src/taxonomy.cpp

Autores originales: Benjamin Callahan, Paul McMurdie y Susan Holmes.
Licencia original LGPL-2, conservada en `COPYING.LIB`. Las modificaciones
del proyecto en `taxonomy_seeded.cpp` se distribuyen bajo la misma licencia.
Fecha de modificación: 2026-10-31. No es una nueva versión oficial de DADA2.

El adaptador `r_seeded_query_mt19937_v1` cambia únicamente la fuente de azar
del desempate. Después del vector bootstrap original, el hilo principal
obtiene una uniforme de R por consulta y la transforma en una semilla uint32
mediante `floor(u * 2^32)`. Cada consulta utiliza su propio `std::mt19937`
para todos sus desempates: secuencia, reverso complementario y 100 réplicas.
La selección uniforme por reservorio entre puntuaciones empatadas se conserva.
No se cambia el cálculo de probabilidades, palabras de ocho bases, orientación,
número de réplicas ni el muestreo bootstrap original.

Se sustituyen además la cabecera interna `dada.h` por cabeceras estándar y
se añaden las anotaciones necesarias para compilar con Rcpp/RcppParallel.
`seeded_classifier.R` clona el entorno de la función R y sustituye únicamente
su entrada nativa; el paquete instalado permanece intacto. Se rechaza otra
versión de DADA2 y se registran ambas fuentes, el cargador y sus hashes.
Cada invocación recompila la biblioteca para evitar reutilizar un binario de
otra configuración local; la caché regenerable permanece fuera de Git.

Contrato: mismos datos **en el mismo orden**, referencia en el mismo orden,
semilla, parámetros y entorno registrado producen tablas iguales al repetir,
incluidos cambios de uno a dos hilos. No se promete invariancia al reordenar
consultas/referencias, otra plataforma, compilador o versión de biblioteca.
Los soportes pueden diferir de ejecuciones históricas cuyo desempate no estaba
controlado. Se conservan antes de reclasificar y se comparan sin seleccionar
la corrida por sus resultados. Las pruebas ambiguas son datos sintéticos,
no evidencia biológica.
