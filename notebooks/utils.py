import operator
import math
import random
import numpy as np
import logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, Polygon
from geopy.distance import great_circle

from simanneal import Annealer

import cartoframes
from cartoframes import *
from cartoframes.viz import *




def location(data: pd.DataFrame, 
             id_col: str, 
             geometry_col: str) -> dict :
    """
    Extract Location from dataframe. Output a dict as {id: (lng, lat), ...}
    """
    loc = {}
    for row in data.iterrows():
        loc_id = row[1][id_col]
        x, y = row[1][geometry_col].x, row[1][geometry_col].y
        loc[loc_id] = loc.get(loc_id, (x,y))
    logging.info(f'[Done] Transform DataFrame To Location Dict)')
    return loc





