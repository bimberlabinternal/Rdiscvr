# R DISCVR
This in an R package that provides helper functions to query data in LabKey/DISCVR modules. This is mainly for internal use by the Bimber Lab.

### <a name="installation">Installation</a>

```{r}
# Install requirements.  Other dependencies will be downloaded automatically
install.packages(c('devtools', 'BiocManager', 'remotes'), dependencies = TRUE, ask = FALSE)

# Updating your Rprofile (i.e. ~/.Rprofile), with the following line will ensure install.packages() pulls from Bioconductor repos:
local({options(repos = BiocManager::repositories())})

# Install latest version:
devtools::install_github(repo = 'bimberlabinternal/Rdiscvr', dependencies = TRUE, upgrade = 'always')
```

### <a name="usage">Usage</a>

Once installed, you will need to set the default baseUrl and defaultFolder:

```{r}
library(Rdiscvr)
    
Rdiscvr::SetLabKeyDefaults(baseUrl = 'https://myserver.com', defaultFolder = 'Labs/Bimber')
```