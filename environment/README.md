# Entorno bioinformático

El entorno completo se ejecuta en Ubuntu/WSL2 porque Bioconda distribuye este
stack para Linux y macOS, no para Windows nativo. El archivo
`environment.yml` fija las herramientas científicas y utiliza únicamente
`conda-forge` y `bioconda` con prioridad estricta.

`r-rcppparallel=5.1.9` está fijado porque el binario Bioconda de DADA2 1.34.0
requiere la ABI TBB que esa versión proporciona. Versiones posteriores pueden
resolver el entorno pero fallan al cargar `dada2.so`; el verificador impide que
ese estado se acepte silenciosamente.

`r-cvxr=1.0_15` está fijado porque ANCOMBC 2.8.0 importa `solve()` desde CVXR.
CVXR 1.9 resolvía las dependencias, pero retiró esa exportación y ANCOMBC no
podía cargarse. Ambas restricciones deben reevaluarse junto con cualquier
actualización deliberada de DADA2 o ANCOMBC.

`r-lme4=1.1_35.5` y `r-matrix=1.6_5` conservan una ABI binaria coincidente para
la generación de ANCOMBC declarada. `r-reformulas` se incluye explícitamente
porque `lmerTest` lo importa aunque el paquete resuelto no lo había instalado.

## Crear o actualizar

Desde PowerShell, abrir WSL y ejecutar en la raíz Linux del repositorio:

```bash
micromamba create --name cacao-microbiome \
  --file environment/environment.yml \
  --channel-priority strict
```

Para reconstruir exactamente el conjunto de binarios ya validado en `linux-64`:

```bash
micromamba create --name cacao-microbiome \
  --file environment/conda-linux-64.lock
```

El lock explícito es la vía recomendada para reproducir este hito. La creación
desde `environment.yml` resuelve de nuevo los paquetes y se reserva para una
actualización deliberada seguida por la validación completa.

Para actualizar deliberadamente un entorno existente, incluida la incorporación
de nuevas dependencias declaradas:

```bash
MAMBA_CHANNEL_PRIORITY=strict micromamba env update \
  --name cacao-microbiome \
  --file environment/environment.yml
```

## Usar sin activación interactiva

```bash
micromamba run --name cacao-microbiome fastqc --version
micromamba run --name cacao-microbiome snakemake --version
micromamba run --name cacao-microbiome Rscript -e 'packageVersion("dada2")'
```

La comprobación reproducible completa se ejecuta dentro del entorno:

```bash
micromamba run --name cacao-microbiome \
  python scripts/environment/verify_environment.py
```

El comando carga los paquetes R requeridos —incluidas las dependencias binarias
críticas— y actualiza
`environment/software_versions.tsv` de forma atómica.

La ubicación física del entorno es una configuración local y no debe
versionarse. Las versiones resueltas y las pruebas de carga se registran en los
artefactos pequeños de `environment/`.
