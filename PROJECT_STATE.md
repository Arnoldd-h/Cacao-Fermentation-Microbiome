# Estado del proyecto

Última actualización: 2026-08-19

## Fase actual

Inicialización del repositorio y definición de las normas de reproducibilidad.

## Completado

- Repositorio remoto renombrado a `Cacao-Fermentation-Microbiome`.
- Rama estable `main` clonada y conectada con `origin`.
- Política de commits, trazabilidad y manejo de datos incorporada al repositorio.
- Exclusiones iniciales para datos ómicos pesados, credenciales, entornos y cachés.

## En progreso

- Definición del alcance científico y de la estructura reproducible del estudio.

## Bloqueado

- Sin bloqueos registrados.

## Decisiones importantes

- `main` debe conservar un estado razonablemente funcional y reproducible.
- Los datos crudos no se versionarán; se reconstruirán desde accessions,
  manifests y scripts.
- Los cambios metodológicos deberán quedar documentados y validados.
- Los commits locales se crearán por unidades lógicas; los pushes requerirán
  autorización explícita.

## Próximas tareas

1. Definir la pregunta de investigación y los objetivos.
2. Documentar criterios de inclusión y exclusión de estudios y muestras.
3. Diseñar la estructura de código, configuración, metadata, tests y resultados.
4. Implementar el inventario inicial de datasets públicos.

## Limitaciones conocidas

- Todavía no existen datasets integrados, pipelines analíticos ni tests.
- Los criterios científicos de inclusión y armonización aún no están definidos.
