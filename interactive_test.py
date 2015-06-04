import karta
import karta_map as km
import matplotlib.pyplot as plt

ax = plt.axes()

ll = (1404428.2170640533, -3011806.030852761)
lr = (2722173.529195537, -1906086.4254576433)
ur = (1975578.9630772884, -1383315.2822752385)
ul = (1019243.9280694319, -2185775.6574322423)

ax.set_xlim(ll[0], lr[0])
ax.set_ylim(ll[1], ur[1])

km.add_graticule(ax, range(-20, 20, 5), range(54, 72, 4), ax_crs=karta.crs.NSIDCNorth)
km.label_ticks(ax, range(-20, 20, 5), range(54, 72, 4), ax_crs=karta.crs.NSIDCNorth)

plt.show()
