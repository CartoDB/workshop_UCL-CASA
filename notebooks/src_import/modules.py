#%matplotlib inline
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from matplotlib.ticker import FormatStrFormatter
import matplotlib.colors as colors
import matplotlib.cm as cm
import matplotlib as mpl
import seaborn as sns
from matplotlib.patches import Ellipse

import importlib
import os
import pathlib
from decimal import Decimal

import pandas as pd
from dask import dataframe as dd
from netCDF4 import Dataset
import numpy as np
from skgstat import Variogram

import cartoframes
from cartoframes import *
from cartoframes.viz import * 
from cartoframes.viz.helpers import * 
from cartoframes.auth import set_default_context 
from cartoframes.viz import Map, Layer, basemaps
from carto.auth import APIKeyAuthClient
from carto.sql import SQLClient
from carto.exceptions import CartoException
from cartoframes.viz.helpers import color_continuous_layer

import geopandas as gpd
from geopandas import GeoDataFrame
import fiona
import shapely
import shapely.wkt
from shapely.geometry import shape, mapping, Point, Polygon, box
from skgstat import Variogram
from pyproj import Proj
import utm
from datetime import datetime

import pysal
from pointpats import PointPattern, PoissonPointProcess, as_window, G, F, J, K, L, Genv, Fenv, Jenv, Kenv, Lenv
from pointpats.centrography import hull, mbr, mean_center, weighted_mean_center, manhattan_median, std_distance,euclidean_median,ellipse
import utm
from pointpats import PointPattern, as_window
from pointpats import PoissonPointProcess as csr
import pointpats.quadrat_statistics as qs

import ipywidgets as widgets
from IPython.display import display, HTML, clear_output,Image

import rpy2.rinterface

import warnings
warnings.filterwarnings('ignore')
# np.set_printoptions(threshold=np.nan)
pd.options.mode.chained_assignment = None
