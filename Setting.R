# pip install --trusted-host pypi.org --trusted-host files.pythonhosted.org os


# �w�� reticulate �]
install.packages("reticulate")
# ���J reticulate �]
library(reticulate)
Sys.which("python")
Sys.getenv("RETICULATE_PYTHON")
# �ˬd�z���t�άO�_�w�˹L Python
py_available(TRUE)
py_available()
conda_install()
# use_condaenv("C:/Users/12476/AppData/Local/Continuum/anaconda3/env")
use_condaenv("env_ret")
use_python("C:\\Users\\12476\\AppData\\Local\\Continuum\\anaconda3\\python.exe", required = T)
# �ˬd���|
py_discover_config()
py_config()


# python �w�� module
# �s�W�����ܼ� C:\Users\12476\AppData\Roaming\Python\Python37\Scripts
# pip install --trusted-host pypi.org --trusted-host files.pythonhosted.org os


# Conda �t�C
conda_list()
conda_create('r-reticulate', packages = "python=3.7")
conda_version()
conda_create("emajor")
conda_info()
py_install("pandas")


# Terminal
# python -m pip install --trusted-host pypi.org --trusted-host files.pythonhosted.org --upgrade pip ����s pip����
# python -m pip --version �ˬd pip ����
# python -m pip install beautifulsoup4
# python -m pip list�C�X�w�U�����M��C


# os.listdir()
# os.environ['R_HOME']
# os.environ['R_USER']