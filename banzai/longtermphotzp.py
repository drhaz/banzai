import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import dateutil.parser



data = np.genfromtxt("photzp.db", unpack=True, dtype=None,\
                  converters={ 1: lambda x: dateutil.parser.parse(x)}, names = ['name','dateobs','airmass','filter','zp'])


plt.plot (data['dateobs'], data['zp'], ".")
plt.gcf().autofmt_xdate()
plt.xlabel ("DATE-OBS")
plt.ylabel ("Photometric Zeropint g")
plt.savefig ('photzptrend.png')


plt.figure()

plt.plot (data['airmass'], data['zp'], ".")
plt.xlabel ("Airmass")
plt.ylabel ("Photomertic Zeropint g")
plt.savefig ('airmasstrend.png')