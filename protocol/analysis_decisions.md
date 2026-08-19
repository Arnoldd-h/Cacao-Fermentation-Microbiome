# Decisiones de análisis

| Fecha | Decisión | Justificación | Alternativas consideradas | Impacto esperado | Estado |
|---|---|---|---|---|---|
| 2026-08-19 | Limitar la primera fase a amplicones bacterianos 16S rRNA. | Evita mezclar tecnologías y permite estabilizar un pipeline defendible. | Integrar ITS o shotgun desde el inicio. | Menor alcance inicial, mayor comparabilidad. | Adoptada |
| 2026-08-19 | Procesar ASVs por estudio y armonizar después a género o familia. | Primers, regiones y protocolos diferentes no producen ASVs directamente comparables. | Concatenar tablas ASV cross-study. | Reduce falsa equivalencia y conserva efectos dentro de estudio. | Adoptada |
| 2026-08-19 | Mantener tiempo absoluto y relativo; derivar etapas desde configuración. | Las fermentaciones tienen duraciones distintas y los límites deben someterse a sensibilidad. | Hardcodear `early/mid/late` en scripts. | Facilita modelos continuos y análisis de sensibilidad. | Adoptada |
| 2026-08-19 | Usar PRJNA492720 como candidato piloto tras verificar su subconjunto publicado de 94 runs 16S. | Tiene diseño longitudinal espontáneo, V4, paired-end y metadata temporal recuperable. | Empezar con todos los proyectos o con WGS. | Permite una vertical slice contenida antes de descargas masivas. | Adoptada |
| 2026-08-19 | Mantener PRJNA865318 y PRJNA1104253 como pendientes. | El primero tiene tiempo incompleto/ambiguo; el segundo mezcla diseños y su subconjunto 16S verificable es controlado in vitro. | Incluirlos automáticamente por contener `AMPLICON`. | Evita ampliar el análisis primario sin evidencia suficiente. | Adoptada |
| 2026-08-19 | Excluir PRJNA552479, PRJNA1257864 y PRJNA1264670 de Fase I. | ENA reporta WGS, no amplicones bacterianos 16S. | Procesarlos como si fueran 16S o ampliar ahora a shotgun. | Mantiene el alcance y deja trazabilidad para fases futuras. | Adoptada |

## Pendiente

- Versión de SILVA y método definitivo de clasificación.
- Umbrales de QC derivados del piloto.
- Tratamiento de ceros y pseudoconteo para CLR.
- Método definitivo de abundancia diferencial según el diseño disponible.
- Selección de licencia del repositorio.
