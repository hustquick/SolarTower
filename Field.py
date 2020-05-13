import numpy as np
from Mirror import Mirror
from matplotlib import pyplot
from shapely.geometry import Polygon, MultiPolygon
from shapely.ops import unary_union
from descartes import PolygonPatch
import ephem
import datetime
import geocoder


# 将向量单位化
def uni_vector(v):
    length = np.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    return v/length


# 求两向量的角平分线向量
def angular_bisector_vector(v_1, v_2):
    v1_uni = uni_vector(v_1)
    v2_uni = uni_vector(v_2)
    return uni_vector(v1_uni + v2_uni)


# 将向量的高度角和方位角形式改写成坐标形式
def angle2vector(alt, azi):
    angle1 = np.pi/2 - np.deg2rad(alt)
    x = np.sin(angle1)*np.sin(np.deg2rad(azi))
    y = np.sin(angle1)*np.cos(np.deg2rad(azi))
    z = np.cos(angle1)
    return np.array([x, y, z])


# 求两个向量之间的夹角的余弦值
def cosine(v_1, v_2):
    uni_v_1 = uni_vector(v_1)
    uni_v_2 = uni_vector(v_2)
    return uni_v_1.dot(uni_v_2)


# 根据向量旋转前的位置和旋转后的位置，找出其旋转轴和旋转角（右手方向）
def calculate_rotation(vector_before, vector_after):
    rotation_axis_ = np.cross(vector_before, vector_after)
    rot_axis = uni_vector(rotation_axis_)
    rot_angle = np.arccos(np.dot(vector_before, vector_after))
    return rot_axis, rot_angle


# 求解绕旋转轴unit_rotation_axis旋转rotation_angle角度时的旋转矩阵
def axis_rotation_matrix(unit_rotation_axis, rotation_angle):
    x = unit_rotation_axis[0]
    y = unit_rotation_axis[1]
    z = unit_rotation_axis[2]
    A_1 = np.array([
        [x*x, x*y, x*z],
        [y*x, y*y, y*z],
        [z*x, z*y, z*z],
    ])
    A_2 = np.array([
        [0, -z, y],
        [z, 0, -x],
        [-y, x, 0],
    ])
    I3 = np.eye(3)
    T = A_1 + (I3 - A_1) * np.cos(rotation_angle) + A_2 * np.sin(rotation_angle)
    return T


# 求解一个点在另一个坐标系中的坐标值。该坐标系给定了z轴的方向，x轴位于原坐标系x0y0平面上。
def coordinate_rotation_matrix(point, vector_z):
    """
    vector_x is defined in the xy plane
    """
    z1, z2, z3 = vector_z
    vector_x = np.array([-z2, z1, 0])/np.sqrt(z1*z1 + z2*z2)
    vector_y = np.cross(vector_z, vector_x)
    x1, x2, x3 = vector_x
    y1, y2, y3 = vector_y
    R = np.array([
        [x1, y1, z1],
        [x2, y2, z2],
        [x3, y3, z3],
                 ])
    p1, p2, p3 = point
    p4 = np.array([p1, p2, p3]).T
    p4_p = np.dot(R, p4)
    point_p = p4_p.ravel()
    return point_p


def alpha(distance):
    if distance < 1000:
        return 0.99321-0.0001176*distance + 1.97e-8*distance**2
    else:
        return np.exp(-0.0001106*distance)


# 获取当前IP地址的经纬度
def get_lalo():
    g = geocoder.ip('me')
    latlng = g.latlng
    print("\n根据你的IP地址")
    # print("你所在城市为:\n" + g.city)
    print("你的经度为：\t%8.4f \n你的纬度为：\t%8.4f" % (latlng[1], latlng[0]))
    return latlng


mirror_typical = Mirror()   # 取一面典型的镜子，获取其参数值。
mirror_length = mirror_typical.length
mirror_center_height = mirror_typical.center_height
tower_height = 115
min_distance = 11.6  # Should be larger than mirror_length * sqrt(2)
base_radius = 115
number_of_rings = 16

# 1代表各镜面坐标系，2代表塔坐标系（世界坐标系，该坐标系中取正东为x轴正方向，正北为y轴正方向，竖直向上为z轴正方向），
# 3代表光线坐标系
v1 = np.array([0, 0, 1])     # 镜面法线向量在镜面坐标系中的坐标值
# 计算太阳高度角
lalo = [37.22, 97.23]   # 手动指定纬度和经度，也可以注释该行，反注释下行，自动获取当前位置的经纬度进行计算
# lalo = get_lalo()     # 计算当前经纬度
observer = ephem.Observer()
observer.lat, observer.lon = str(lalo[0]), str(lalo[1])

observer.date = ephem.Date('2020/03/20 06:00:00.00')  # 手动指定时间(UTC时间)，也可以注释该行，反注释下面两行，自动获取当前时间进行计算
# gatech.date = datetime.datetime.utcnow()    # 计算当前时间
# print("\n当前时间为:\n%s" % ephem.localtime(gatech.date))

sun = ephem.Sun()
sun.compute(observer)

print("太阳高度角: %s \n太阳方位角: %s" % (sun.alt, sun.az))
altitude = np.rad2deg(sun.alt)  # 也可以注释该两行，反注释下面两行，手动输入需要计算的高度角和方位角
azimuth = np.rad2deg(sun.az)
# altitude = 52.78
# azimuth = 180
v3 = angle2vector(altitude, azimuth)    # 光线向量在世界坐标系中的坐标值，指向太阳！

radius = np.zeros(number_of_rings)      # 用来记录各圈的半径
number = np.zeros(number_of_rings, np.int)  # 用来记录各圈的镜子数目
radius[0] = base_radius
for i in range(number_of_rings - 1):
    number[i] = np.floor(2*np.pi*radius[i]/min_distance)
    distance_to_increase = np.ceil(mirror_length * radius[i] / 105)  # Why?
    radius[i+1] = radius[i] + distance_to_increase
number[number_of_rings-1] = np.floor(2*np.pi*radius[number_of_rings-1]/min_distance)

# 定义镜子列表
mirror_list = [[Mirror() for j in range(number[i])] for i in range(number_of_rings)]
number_of_mirrors = sum([len(mirror_list[i]) for i in range(len(mirror_list))])

points2 = []    # 用于记录各镜子顶点旋转后在世界坐标系中的坐标值
points3 = []    # 用于记录各镜子顶点旋转后在光线坐标系中的坐标值
for i in range(number_of_rings):
    points2.append(np.zeros((number[i], 4, 3)))
    points3.append(np.zeros((number[i], 4, 3)))
    for j in range(number[i]):
        # 计算镜子中心坐标
        delta_theta = 2 * j * np.pi / number[i] + np.pi/2  # 记录镜子与塔之间连线与正东方向的夹角，逆时针为正
        mirror_list[i][j].center = np.array([radius[i] * np.cos(delta_theta),
                                            radius[i] * np.sin(delta_theta),
                                            mirror_center_height])

        # 计算镜子法线向量
        v_m2t = np.array([0, 0, tower_height]) - mirror_list[i][j].center   # 镜子到塔顶的向量
        mirror_list[i][j].orientation = angular_bisector_vector(v3, v_m2t)

        # 计算镜子余弦损失
        v2 = mirror_list[i][j].orientation  # 旋转之后，镜子的法线在世界坐标系中的坐标值
        mirror_list[i][j].cosine_efficient = cosine(v2, v3)
        # 计算镜子旋转轴及旋转角
        rotation_axis, rotation_angle = calculate_rotation(v1, v2)
        # 计算绕轴旋转的旋转矩阵
        T = axis_rotation_matrix(rotation_axis, rotation_angle)
        # 计算旋转之后的镜子各顶点坐标
        for k in range(4):
            # 各顶点在世界坐标系中的坐标值
            points2[i][j][k] = np.dot(mirror_list[i][j].points[k].co, T) + mirror_list[i][j].center
            # 各顶点在光线坐标系中的坐标值
            points3[i][j][k] = coordinate_rotation_matrix(points2[i][j][k], v3)

# 计算各镜面在光线方向的投影面积及遮挡之后的总面积
GM = (np.sqrt(5) - 1.0) / 2.0
W = 8.0
H = W * GM
SIZE = (W, H)
BLUE = '#6699cc'
GRAY = '#999999'

polygons = []
for i in range(number_of_rings):
    for j in range(number[i]):
        points_to_sunshine = [
            points3[i][j][0][0:2],
            points3[i][j][1][0:2],
            points3[i][j][2][0:2],
            points3[i][j][3][0:2],
        ]
        polygons.append(Polygon(points_to_sunshine))

fig = pyplot.figure(1, figsize=SIZE, dpi=300)

ax = fig.add_subplot(121)

for ob in polygons:
    p = PolygonPatch(ob, fc=GRAY, ec=GRAY, alpha=0.5, zorder=1)
    ax.add_patch(p)

ax.set_title('a) Before overlap')
ax.axis('equal')
# set_limits(ax, -2, 6, -2, 2)
ax = fig.add_subplot(122)

u = unary_union(polygons)
patch2b = PolygonPatch(u, fc=BLUE, ec=BLUE, alpha=0.5, zorder=2)
ax.add_patch(patch2b)

ax.set_title('b) After overlap')
ax.axis('equal')
# set_limits(ax, -2, 6, -2, 2)

pyplot.savefig('picture.pdf')
pyplot.show()

m = MultiPolygon(polygons)
print("迎着太阳的总面积为：{:.2f}".format(m.area))
print("被遮挡之后的总面积为：{:.2f}".format(u.area))