install.packages("caret")
install.packages("mt")
install.packages('ROCR')
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("preprocessCore")

install.packages('zoo')
install.packages("GSE")
install.packages("remotes")
remotes::install_github("bioinf-jku/platt")

#如果rtools44没有配置成功, 先在conda环境中的R里边配置变量，然后重启就可以了
#writeLines('PATH="C:\\rtools44\\usr\\bin;C:\\rtools44\\mingw64\\bin;${PATH}"', "~/.Renviron")
