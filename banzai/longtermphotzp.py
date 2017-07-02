import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime
import dateutil.parser



data = np.genfromtxt("photzp.db", unpack=True, dtype=None,\
                  converters={ 1: lambda x: dateutil.parser.parse(x)}, names = ['name','dateobs','airmass','filter','zp'])


meanzp = np.nanmedian(data['zp'])
zp_air = data['zp'] + 0.2 * data['airmass'] - 0.2

plt.plot (data['dateobs'], data['zp'], ".", c="grey")
plt.plot (data['dateobs'], zp_air, ".", c="blue")

plt.ylim([meanzp-0.5,meanzp+0.5])
plt.gcf().autofmt_xdate()
plt.xlabel ("DATE-OBS")
plt.ylabel ("Photometric Zeropint g")
plt.savefig ('photzptrend.png')


plt.figure()

plt.plot (data['airmass'], data['zp'], ".", c="grey")
plt.plot (data['airmass'], zp_air, ".", c="blue")

plt.xlabel ("Airmass")
plt.ylabel ("Photomertic Zeropint g")
plt.ylim([meanzp-0.5,meanzp+0.5])

plt.savefig ('airmasstrend.png')