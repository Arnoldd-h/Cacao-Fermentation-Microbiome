# Resúmenes descriptivos por familia y género

## Política registrada

`61da448` fija `config/taxonomic_reporting.yaml` y la decisión correspondiente
antes de generar estas tablas. Se consumen exclusivamente las tablas derivadas
de la separación bacteriana, cuya validación independiente debe pasar de nuevo.
Los perfiles `pilot` y `full_study` esperan seis y 93 muestras respectivamente.

Cada ASV se agrega únicamente con ASVs de la misma muestra que comparten todos
los rangos originales hasta familia o género. Dos etiquetas de género iguales
con familias distintas permanecen separadas. Las asignaciones vacías conservan
sus ancestros conocidos; las de Kingdom desconocido también se retienen.
No se completa una etiqueta a partir de otro rango ni se alteran llamadas SILVA.

La abundancia relativa es `lecturas del linaje / todas las lecturas retenidas
de la muestra`. El denominador incluye lecturas sin asignación suficiente.
Las filas de cero se omiten de la representación larga, sin filtrar ASVs,
lecturas ni muestras. No se añaden pseudoconteos. La cobertura separa etiquetas
vacías y no vacías; una etiqueta no vacía no garantiza identificación de especie
ni corresponde necesariamente a un género cultivado.

## Productos y validación

- `sample_composition.tsv`: estudio, muestra, run, lote, horas, tiempo relativo,
  etapa, estrato, rango, linaje, conteo, total retenido y abundancia relativa.
- `sample_label_coverage.tsv`: lecturas con/sin etiqueta en cada rango y su
  fracción respecto al mismo total retenido.
- `summary.json`, configuración, hashes de entradas/salidas, procedencia,
  marcador `SUCCESS` y `validation.json`.

El validador recorre los conteos fuente, comprueba la pertenencia exacta de
muestras/ASVs, reconcilia cada linaje y verifica conteos, metadata, denominadores,
cobertura y hashes. La tolerancia numérica procede de la configuración de
diversidad existente. Las pruebas sintéticas cubren homónimos con ancestros
distintos, etiquetas vacías, Kingdom desconocido, ceros, estudios mezclados,
identidades duplicadas y alteraciones de salidas. No representan datos observados.

## Ejecución

Después de completar y validar la separación del perfil elegido:

```bash
python3 scripts/environment/run_in_environment.py python scripts/reporting/summarize_composition.py --analysis full_study
python3 scripts/environment/run_in_environment.py python scripts/reporting/summarize_composition.py --analysis full_study --validate-only
```

Para el piloto se utiliza `--analysis pilot`. No hay nuevas dependencias.
Las tablas conservan muestras y estratos individuales. No agregan lotes o
tiempos, no estiman asociaciones o significación y no definen un core. Son
entradas descriptivas para revisar la composición antes de registrar y ejecutar
los modelos longitudinales del protocolo.

## Primera ejecución observada

El perfil piloto ejecutado sobre `ad4c3ca` conserva seis muestras, 88 ASVs y
123.693 lecturas en 95 filas de composición. La validación de 15 entradas y
cinco artefactos pasa tanto al generar como al comprobar las salidas existentes.
El perfil completo contiene 93 muestras, 764 ASVs y 2.096.980 lecturas en
1.477 filas de composición. Se conservan los cuatro lotes y sus estratos;
15 entradas y cinco artefactos pasan validación. El hash de la ejecución
vigente figura en su procedencia, sin equiparar los ASVs del perfil completo
con los identificadores del piloto.
