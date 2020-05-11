from matplotlib import pyplot
from shapely.geometry import Polygon, MultiPolygon
from shapely.ops import unary_union
from descartes import PolygonPatch
from math import sqrt

GM = (sqrt(5) - 1.0) / 2.0
W = 8.0
H = W * GM
SIZE = (W, H)

BLUE = '#6699cc'
GRAY = '#999999'

mirror_list = [
                [[0, 0], [4, 0], [4, 4], [0, 4]],
                [[6, 0], [10, 0], [10, 4], [6, 4]],
                [[3, 3], [7, 3], [7, 7], [3, 7]],
                [[7, 1], [3, 1], [3, -3], [7, -3]],
               ]

polygons = []
for mirror in mirror_list:
    polygons.append(Polygon(mirror))

fig = pyplot.figure(1, figsize=SIZE, dpi=300)
ax = fig.add_subplot(121)

for ob in polygons:
    p = PolygonPatch(ob, fc=GRAY, ec=GRAY, alpha=0.5, zorder=1)
    ax.add_patch(p)

ax.set_title('a) Before overlap')

ax.axis("equal")
ax = fig.add_subplot(122)

u = unary_union(polygons)
patch2b = PolygonPatch(u, fc=BLUE, ec=BLUE, alpha=0.5, zorder=2)
ax.add_patch(patch2b)

ax.set_title('b) After overlap')

ax.axis("equal")

pyplot.savefig('picture.pdf')
pyplot.show()

m = MultiPolygon(polygons)
print("The area before overlap：{:.2f}".format(m.area))
print("The area after overlap：{:.2f}".format(u.area))
