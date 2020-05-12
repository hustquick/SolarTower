import numpy as np
from Mirror import Mirror
from matplotlib import pyplot
from shapely.geometry import Polygon, MultiPolygon
from shapely.ops import unary_union
from descartes import PolygonPatch


def uni_vector(v):
    length = np.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    return v/length


def angular_bisector_vector(v1, v2):
    v1_uni = uni_vector(v1)
    v2_uni = uni_vector(v2)
    return uni_vector(v1_uni + v2_uni)


def angle2vector(alt, azi):
    angle1 = np.pi/2 - alt
    x = np.sin(angle1)*np.sin(azi)
    y = np.sin(angle1)*np.cos(azi)
    z = np.cos(angle1)
    return np.array([x, y, z])


def vector2angle(a, b, c):
    d = np.sqrt(a**2 + b**2)
    alt = np.arctan2(c, d)
    azi = np.pi/2 - np.arctan2(b, a)
    return alt, azi


def cosine(v1, v2):
    """
    v1 and v2 are unit vectors
    :param v1:
    :param v2:
    :return:
    """
    return v1.dot(v2)


def calculate_rotation(vector_before, vector_after):
    rotation_axis_ = np.cross(vector_before, vector_after)
    rotation_axis = uni_vector(rotation_axis_)
    rotation_angle = np.arccos(np.dot(vector_before, vector_after))
    return rotation_axis, rotation_angle


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
    I = np.eye(3)
    T = A_1 + (I - A_1) * np.cos(rotation_angle) + A_2 * np.sin(rotation_angle)
    return T


def coordinate_rotation_matrix(point, vector_z):
    """
    vector_x is defined in the xy plane
    :param point:
    :param vector_z:
    :return:
    """
    z1, z2, z3 = vector_z
    vector_x = np.array([z2, -z1, 0])/np.sqrt(z1*z1 + z2*z2)
    vector_y = np.cross(vector_z, vector_x)
    x1, x2, x3 = vector_x
    y1, y2, y3 = vector_y
    R = np.array([
        [x1, y1, z1, 0],
        [x2, y2, z2, 0],
        [x3, y3, z3, 0],
        [0, 0, 0, 1],
                 ])
    p1, p2, p3 = point
    p4 = np.array([p1, p2, p3, 1])
    p4_p = np.dot(p4, R)
    point_p = p4_p[0:3]
    return point_p


def alpha(distance):
    if distance < 1000:
        return 0.99321-0.0001176*distance + 1.97e-8*distance**2
    else:
        return np.exp(-0.0001106*distance)


mirror_typical = Mirror()
tower_height = 115
mirror_length = mirror_typical.length
mirror_center_height = mirror_typical.center_height
min_distance = 11.6  # Should be larger than mirror.length * sqrt(2)
base_radius = 115
number_of_rings = 16

# 1代表各镜面坐标系，2代表塔坐标系（世界坐标系），3代表光线坐标系
v1 = np.array([0, 0, 1])     # 镜面法线向量在镜面坐标系中的坐标值
# zenith = np.deg2rad(37.22)
altitude = np.deg2rad(52.78)
azimuth = np.deg2rad(0)
v3 = angle2vector(altitude, azimuth)    # 光线向量在世界坐标系中的坐标值


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
        delta_theta = 2 * j * np.pi / number[i]  # 记录镜子与塔之间连线与正北方向的夹角，逆时针为正
        mirror_list[i][j].center = np.array([radius[i] * np.cos(delta_theta-np.pi/2),
                                            radius[i] * np.sin(delta_theta-np.pi/2),
                                            mirror_center_height])

        # 计算镜子法线向量
        v_m2t = np.array([0, 0, tower_height]) - mirror_list[i][j].center
        mirror_list[i][j].orientation = angular_bisector_vector(v3, v_m2t)

        # 计算镜子余弦损失
        mirror_normal_angle = vector2angle(mirror_list[i][j].orientation[0],
                                           mirror_list[i][j].orientation[1],
                                           mirror_list[i][j].orientation[2]
                                           )
        mirror_list[i][j].cosine_efficient = cosine(mirror_list[i][j].orientation, v3)

        v2 = mirror_list[i][j].orientation  # 旋转之后，镜子的法线在世界坐标系中的坐标值
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