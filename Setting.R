# pip install --trusted-host pypi.org --trusted-host files.pythonhosted.org os


# 安裝 reticulate 包
install.packages("reticulate")
# 載入 reticulate 包
library(reticulate)
Sys.which("python")
Sys.getenv("RETICULATE_PYTHON")
# 檢查您的系統是否安裝過 Python
py_available(TRUE)
py_available()
conda_install()
# use_condaenv("C:/Users/12476/AppData/Local/Continuum/anaconda3/env")
use_condaenv("env_ret")
use_python("C:\\Users\\12476\\AppData\\Local\\Continuum\\anaconda3\\python.exe", required = T)
# 檢查路徑
py_discover_config()
py_config()


# python 安裝 module
# 新增環境變數 C:\Users\12476\AppData\Roaming\Python\Python37\Scripts
# pip install --trusted-host pypi.org --trusted-host files.pythonhosted.org os


# Conda 系列
conda_list()
conda_create('r-reticulate', packages = "python=3.7")
conda_version()
conda_create("emajor")
conda_info()
py_install("pandas")


# Terminal
# python -m pip install --trusted-host pypi.org --trusted-host files.pythonhosted.org --upgrade pip 先更新 pip版本
# python -m pip --version 檢查 pip 版本
# python -m pip install beautifulsoup4
# python -m pip list列出已下載的套件。


# os.listdir()
# os.environ['R_HOME']
# os.environ['R_USER']