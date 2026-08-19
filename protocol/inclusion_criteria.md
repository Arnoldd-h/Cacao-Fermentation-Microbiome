# Criterios de inclusión y exclusión

## Inclusión primaria

Un estudio puede incluirse en el análisis principal cuando cumple todos los
criterios siguientes:

1. Contiene amplicones bacterianos 16S rRNA identificables sin ambigüedad.
2. Estudia fermentación de granos o masa pulpa-grano de cacao.
3. Corresponde a una fermentación espontánea; los starters se analizan aparte.
4. Tiene FASTQ públicos y accessions de run recuperables.
5. Dispone de al menos tres tiempos u orden temporal inequívoco para sucesión.
6. Permite conservar `study_id`, muestra, run y lote de fermentación.
7. La región 16S y los primers están publicados o pueden verificarse antes del
   procesamiento definitivo.
8. La metadata permite separar muestras de fermentación, controles, otros
   marcadores y otros ambientes.

## Inclusión condicional

Se asigna `pending` cuando el estudio parece pertinente pero falta verificar
alguno de estos puntos:

- paper o protocolo de laboratorio;
- primers/región 16S;
- significado de identificadores de muestra;
- duración o tiempo de fermentación;
- separación entre fermentación espontánea y tratamiento inoculado;
- independencia respecto de otro BioProject o publicación.

Los estudios pendientes no aportan evidencia al análisis principal.

## Exclusión primaria

Se excluyen de la Fase I:

- ITS, shotgun/WGS o datos no 16S sin un subconjunto 16S separable;
- phyllosphere, suelo, fruto intacto u otras matrices no correspondientes a la
  masa de fermentación;
- mock communities, controles negativos y controles técnicos;
- estudios de un solo tiempo sin contraste temporal;
- datos sin reads crudos públicos;
- muestras con tiempo no recuperable tras revisar fuentes primarias;
- fermentaciones exclusivamente inoculadas para el análisis de sucesión
  espontánea principal;
- duplicados del mismo conjunto de runs.

## Criterios posteriores de QC

Los umbrales de lecturas y retención se definirán después del piloto, usando la
distribución real de calidad. No se fija todavía un corte arbitrario. Toda
exclusión posterior deberá indicar `sample_id`, razón, métrica, umbral y decisión
en el registro de exclusiones.

## Resolución de conflictos

Ante discrepancias, tienen prioridad los registros de run/experimento y el
método publicado sobre inferencias desde títulos. Las inferencias necesarias se
marcan explícitamente y no se presentan como metadata depositada.
