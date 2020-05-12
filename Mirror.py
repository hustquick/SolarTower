import numpy as np


class Point:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        self.co = np.array([x, y, z])


class Mirror:
    def __init__(self, length=8, width=8, center_height=6):
        self.length = length
        self.width = width
        self.center_height = center_height
        self.center = None
        self.orientation = None
        self.points = [
            Point(self.length / 2, self.width / 2, 0),
            Point(-self.length / 2, self.width / 2, 0),
            Point(-self.length / 2, -self.width / 2, 0),
            Point(self.length / 2, -self.width / 2, 0),
        ]
        self.cosine_efficient = None

