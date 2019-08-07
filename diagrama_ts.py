# Importe as dependencias
from __future__ import unicode_literals
import matplotlib.pyplot as plt
import iris
import iris.plot as iplt
import iris.quickplot as qplt
import netCDF4
from datetime import datetime, timedelta
import os, sys, string
import numpy as np

_author_   = 'Ueslei Adriano Sutil'
_email_    = 'ueslei@outlook.com'
_created_  = datetime(2017, 03, 20)
_modified_ = datetime(2017, 03, 20)
_version_  = "0.1.0"
_status_   = "Development"

# Abra o arquivo e carregue as variaveis.
nc = '/home/uesleisutil/Documentos/python_scripts/myocean_anita.nc'

theta = iris.load_cubes('/home/uesleisutil/Documentos/python_scripts/myocean_anita.nc')
