intensity = 20
sharpness = 0.5


class Box:
    def __init__(self, x1, x2, y1, y2, z1, z2):
        self.x1 = x1
        self.x2 = x2
        self.y1 = y1
        self.y2 = y2
        self.z1 = z1
        self.z2 = z2

        self.intensity = intensity
        self.sharpness = sharpness

    def point_in_box(self, point: list):
        if len(point) != 3:
            print("No valid point")
            return False
        return (self.x1 <= point[0] <= self.x2) and (self.y1 <= point[1] <= self.y2) and (self.z1 <= point[2] <= self.z2)


# this is a method decoding the boxes for any of the ligand size
def get_boxes(denticity, planar=True):

    box_list = list()

    if denticity == 1:
        pass
    elif denticity == 2:
        box_list.append(Box(-5.0, 5.0, -5.0, 5.0, 4.0, 10.0))
        box_list.append(Box(2.5, 100.0, -5.0, 5.0, -2.5, 5.0))
        box_list.append(Box(0.4, 100.0, -10.0, 10.0, 0.001, 1000.))
        box_list.append(Box(100.0, 101.0, 100.0, 101., 100., 101.))
        box_list.append(Box(-1.0, 1.0, -1.0, 1.0, 0.001, 100.))
    elif denticity == 3:
        box_list.append(Box(0.0, 1.0, -1.0, 1.0, 0.0001, 10.0))
        box_list.append(Box(-1.0, 1.0, -1.0, 1.0, -10.0, -0.2))
        box_list.append(Box(-2.0, 2.0, -2.0, 2.0, 0.8, 10.))
        box_list.append(Box(-6.0, -2.0, -2.0, -1.0, -10.0, 10.0))
        box_list.append(Box(-6.0, -2.0, 1.0, 2.0, -10.0, 10.0))
        box_list.append(Box(-2.0, 2.0, -2.0, 2.0, -10.0, -0.8))
    elif denticity == 5:
        box_list.append(Box(-1.0, 1.0, -1.0, 1.0, 0.3, 2.001))
        box_list.append(Box(100.0, 101.0, 100.0, 101., 100., 101.))
        box_list.append(Box(-1.8, 1.8, -1.8, 1.8, 2.001, 6.001))
        box_list.append(Box(100.0, 101.0, 100.0, 101., 100., 101.))
        box_list.append(Box(-4.0, 4.0, -4.0, 4.0, 6.001, 14.001))
        box_list.append(Box(100.0, 101.0, 100.0, 101., 100., 101.))
    elif denticity == 4 and planar is True:
        box_list.append(Box(-1.0, 1.0, -1.0, 1.0, 0.001, 2.001))
        box_list.append(Box(-1.0, 1.0, -1.0, 1.0, -2.001, -0.001))
        box_list.append(Box(-2.0, 2.0, -2.0, 2.0, -6.001, -2.001))
        box_list.append(Box(-2.0, 2.0, -2.0, 2.0, 2.001, 6.001))
        box_list.append(Box(-4.0, 4.0, -4.0, 4.0, 6.001, 14.001))
        box_list.append(Box(-4.0, 4.0, -4.0, 4.0, -14.001, -6.001))
    elif denticity == 4 and planar is False:
        box_list.append(Box(-10.0, -0.9, -1.8, 1.8, -1.3, 1.3))
        box_list.append(Box(-10.0, -1.8, -5.0, 5.0, -2.0, 2.0))
        box_list.append(Box(-10.0, -3.5, -10.0, 10.0, -10.0, 10.0))
        box_list.append(Box(100.0, 101.0, 100.0, 101., 100., 101.))
        box_list.append(Box(100.0, 101.0, 100.0, 101., 100., 101.))
        box_list.append(Box(100.0, 101.0, 100.0, 101., 100., 101.))
    else:
        pass

    return box_list

