class Point:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z


class Mirror:
    def __init__(self, p1, p2, p3: Point, p4: Point):
        self.p1 = p1
        self.p2 = p2
        self.p3 = p3
        self.p4 = p4
        self.length = 8
        self.center = None
        self.orientation = None

