# Métodos

El método definitivo se completará a partir del protocolo y de las versiones
ejecutadas del workflow. No se describirán análisis que no hayan sido realizados.

## Modelo longitudinal exploratorio dentro de estudio

El plan se fijó en `c97345c`, después del procesamiento y descripción inicial,
antes de ajustar el primer modelo. Se analizaron 764 ASVs bacterianas retenidas
de PRJNA492720, independientemente del piloto. La unidad independiente fue el
lote (cuatro). El principal exigió ambos estratos por lote/hora y utilizó
46 observaciones/92 muestras; la muestra sin pareja se conservó en sensibilidad.

Se transformaron los conteos mediante CLR después de añadir un pseudoconteo
1 a todos los componentes. Dentro de lote/hora se promediaron vectores CLR
con igual peso por estrato. Se ajustó OLS multivariado con intercepto y tiempo
relativo por lote; las pendientes se combinaron con igual peso por lote.
Sensibilidades fijadas: pseudoconteo 0,5, horas/24, estratos disponibles y
cada estrato separado. No se seleccionaron ASVs según abundancia o resultado.

Se enumeraron los 16 cambios de signo de vectores completos de pendiente.
La estadística fue la norma de la pendiente media; p fue la proporción
inclusiva de normas al menos tan grandes. El contraste depende de independencia
y simetría central bajo la nula, no de intercambiabilidad entre horas. Los seis
escenarios formaron una familia FDR BY. Se reportaron ajuste, direcciones,
dispersión de pendientes y leave-one-batch-out; no se realizaron pruebas por ASV.
La resolución mínima p=0,125 se fijó antes de observar los ajustes.

Python y su biblioteca estándar produjeron las estimaciones; R 4.4.3 verificó
independientemente centros, membresía y coeficientes con `stats::lm.fit`,
signos exhaustivos y `stats::p.adjust`. El primer ajuste `7b6d5a6` y su
renovación operativa `6415bfd` conservaron
configuración, versiones, semillas y hashes de 27 entradas/19 artefactos.
Métodos completos, fuentes primarias y límites:
[longitudinal_methods.md](../reports/longitudinal_methods.md).
