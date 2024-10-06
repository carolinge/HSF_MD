import numpy as np
import pandas as pd

input_name="NMR_240526.xlsx"

dataN = pd.read_excel(input_name, skiprows=0, usecols="A:R", index_col=None, sheet_name = "Sheet0", dtype = {"ID":str,"InStore":str})
dataNMR =dataN.interpolate()
dataNMR.to_csv("NMRinter.csv")



