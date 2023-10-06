import os,sys
PROJECT_ROOT = os.path.abspath(os.path.join(
                  os.path.dirname("__file__"), 
                  os.pardir)
)
sys.path.append(PROJECT_ROOT)



import pandas as pd
import numpy as np
from numpy import log as ln
from math import e as e
from math import sin as sin
from math import cos as cos
import openpyxl
import scipy.optimize as opt
#from scipy.misc import derivative
import numdifftools as nd
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots