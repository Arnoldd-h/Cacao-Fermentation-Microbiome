# Política operativa del repositorio

Estas reglas se aplican a cualquier trabajo automatizado o asistido realizado
en este repositorio. El historial Git forma parte del registro reproducible del
estudio.

## Control de versiones

- Antes de trabajar, comprobar si el directorio pertenece a un repositorio Git;
  nunca crear un repositorio anidado.
- Inspeccionar `git status`, el remote y `.gitignore` antes de cambios
  importantes.
- Crear commits locales cuando exista una unidad lógica terminada y verificable;
  no agrupar fases completas ni versionar cada cambio trivial por separado.
- Antes de cada commit, revisar `git status` y `git diff`, ejecutar las pruebas o
  validadores pertinentes y comprobar que no haya credenciales, temporales ni
  datos pesados.
- Preferir commits atómicos y mensajes informativos de Conventional Commits,
  por ejemplo `feat: add SRA metadata inventory pipeline`.
- Mantener `main` razonablemente funcional y reproducible. Usar ramas
  descriptivas para tareas grandes o experimentales, no para cambios triviales.
- No ejecutar `git reset --hard`, `git clean -fd`, force push ni reescrituras
  destructivas del historial compartido.
- No sobrescribir cambios del usuario.

## Estado y trazabilidad científica

- Actualizar `PROJECT_STATE.md` en hitos o cambios importantes: fase actual,
  completado, en progreso, bloqueos, decisiones, tareas siguientes y
  limitaciones.
- Actualizar `CHANGELOG.md` ante nuevas capacidades, cambios metodológicos o de
  arquitectura, nuevas fases, datasets o correcciones científicas relevantes.
- Toda modificación de umbrales, etapas temporales, criterios de exclusión,
  taxonomía, modelos estadísticos o métodos composicionales debe actualizar la
  configuración y `analysis_decisions.md`, ejecutar las validaciones afectadas
  y recibir un commit específico.
- Cuando sea útil, registrar el hash de `git rev-parse HEAD`, la configuración,
  las versiones de software, el dataset y la fecha de ejecución junto a
  resultados científicos importantes.
- Usar tags anotados en hitos consolidados. Secuencia prevista:
  `v0.1.0-dataset-inventory`, `v0.2.0-pilot-qc`, `v0.3.0-pilot-dada2`,
  `v0.4.0-taxonomy`, `v0.5.0-diversity`, `v0.6.0-multistudy`,
  `v0.7.0-meta-analysis`, `v0.8.0-results`, `v0.9.0-manuscript` y `v1.0.0`.

## Datos y resultados

- No versionar FASTQ, BAM, SAM, CRAM, grandes archivos comprimidos, bases de
  datos completas, entornos Conda, cachés, temporales, resultados intermedios
  pesados, credenciales ni claves API.
- Los datos crudos deben reconstruirse mediante accession IDs, manifests y
  scripts de descarga.
- Versionar principalmente código, configuración, metadata, documentación,
  tests, tablas derivadas pequeñas, figuras relevantes y el protocolo
  científico.
- Se permiten resúmenes de QC, tablas estadísticas finales y figuras de
  publicación de tamaño razonable si son trazables al código, configuración y
  versión Git que los generó.

## Commits y repositorio remoto

- Codex está autorizado a crear commits locales periódicos después de revisar y
  validar una unidad lógica de progreso.
- No asumir que cada commit debe enviarse inmediatamente.
- No crear repositorios remotos, cambiar su visibilidad, publicar cambios ni
  hacer push sin autorización explícita del usuario.
- Después de un commit importante, informar hash corto, mensaje, contenido,
  validación y siguiente etapa.

El objetivo es que `git log --oneline --graph` relate de forma comprensible la
evolución desde el inventario y la metadata hasta QC, DADA2, taxonomía,
diversidad, integración, meta-análisis, figuras y manuscrito.
