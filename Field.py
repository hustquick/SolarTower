import numpy as np


def uni_vector(v):
    length = np.sqrt(v[0]**2 + v[1]**2 + v[2]**2)
    return v/length


def angular_bisector_vector(v1, v2):
    v1_uni = uni_vector(v1)
    v2_uni = uni_vector(v2)
    return uni_vector(v1_uni + v2_uni)


def angle2vector(alt, azi):
    theta = np.pi/2 - alt
    x = np.sin(theta)*np.sin(azi)
    y = np.sin(theta)*np.cos(azi)
    z = np.cos(theta)
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
    rotation_axis = np.cross(vector_before, vector_after)
    rotation_axis = uni_vector(rotation_axis)
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


tower_height = 115
mirror_length = 8
mirror_center_height = 6
min_distance = 11.6  # Should be larger than mirror_length * sqrt(2)
base_radius = 115
number_of_rings = 16

zenith = np.deg2rad(37.22)
altitude = np.deg2rad(52.78)
azimuth = np.deg2rad(0)
v_sunshine = angle2vector(altitude, azimuth)

radius = np.zeros(number_of_rings)
number = np.zeros(number_of_rings)
radius[0] = base_radius

for i in range(number_of_rings - 1):
    number[i] = np.floor(2*np.pi*radius[i]/min_distance)
    distance_to_increase = np.ceil(mirror_length * radius[i] / 105)  # Why 105?
    radius[i+1] = radius[i] + distance_to_increase

number_of_mirrors = sum(number)

# mirror = np.zeros(number_of_mirrors)
mirror_center = []
for i in range(len(radius) - 1):
    number_of_mirrors_in_row = int(number[i])
    theta = np.zeros(number_of_mirrors_in_row)
    mirror_center.append(np.zeros([number_of_mirrors_in_row, 3]))
    for j in range(number_of_mirrors_in_row):
        theta[j] = 2 * j * np.pi / number[i]
        mirror_center[i][j] = np.array([radius[i] * np.cos(theta[j]-np.pi/2), radius[i] * np.sin(theta[j]-np.pi/2),
                                        mirror_center_height])

mirror_normal_vector = []
cosine_efficient = []
for i in range(len(mirror_center)):
    mirror_normal_vector.append(np.zeros([len(mirror_center[i]), 3]))
    cosine_efficient.append(np.zeros(len(mirror_center[i])))
    for j in range(len(mirror_center[i])):
        v_m2t = np.array([0, 0, tower_height]) - mirror_center[i][j]
        mirror_normal_vector[i][j] = angular_bisector_vector(v_sunshine, v_m2t)
        mirror_normal_angle = vector2angle(mirror_normal_vector[i][j][0], mirror_normal_vector[i][j][1],
                                           mirror_normal_vector[i][j][2])
        cosine_efficient[i][j] = cosine(mirror_normal_vector[i][j], v_sunshine)



# 1代表各镜面坐标系，2代表塔坐标系（世界坐标系），3代表光线坐标系
points1 = np.array([
    [mirror_length/2, mirror_length/2, 0],
    [-mirror_length/2, mirror_length/2, 0],
    [-mirror_length/2, -mirror_length/2, 0],
    [mirror_length/2, -mirror_length/2, 0],
                           ])
v1_horizontal = np.array([0, 0, 1])
v2 = mirror_normal_vector
rotation_axis, rotation_angle = calculate_rotation(v1_horizontal, v2[0][1])
T = axis_rotation_matrix(rotation_axis, rotation_angle)
points2 = np.dot(points1, T) + mirror_center[0][1]
v3 = v_sunshine
points3 = coordinate_rotation_matrix(points2[0], v3)
