import logging
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
            logging.debug("No valid point")
            return False
        return (self.x1 <= point[0] <= self.x2) and (self.y1 <= point[1] <= self.y2) and (self.z1 <= point[2] <= self.z2)


# this is a method decoding the boxes for any of the ligand size
def get_boxes(denticity, planar: bool = True, input_topology=None, bool_placed_boxes=None, build_options: str = None):
    """
    planar only important for tetradentates
    """

    box_list = list()
    logging.debug("the type of input is: "+str(type(input_topology)))
    logging.debug("The input topology is: "+str(input_topology))
    logging.debug("we are in the box list section of the code")
    if (denticity == 2) and ((input_topology == [2, 2]) or (input_topology == [2, 1, 1] or input_topology == [2, 1, 0]) or (input_topology == [2, 0])) and bool_placed_boxes == False:
        if build_options == 'slab':
            logging.debug("Bidenatate Rotator: Right --> slab chosen")
            box_list.append(Box(-100.0, 100.0, -1.0, 1.0, 0.5, 1.0))        # Top_Plate
            box_list.append(Box(-100.0, 100.0, -1.0, 1.0, -1.0, -0.5))      # Bottom_Plate
            box_list.append(Box(-100.0, -0.1, -1.0, 1.0, -1.0, 1.0))
            box_list.append(Box(-100.0, -0.5, -100, 100, -100, 100))
        elif build_options == 'horseshoe':
            logging.debug("Bidenatate Rotator: Left --> horseshoe chosen")
            box_list.append(Box(-4.0, 4.0, -1.4, 1.4, 0.5, 100))            # Top_Plate
            box_list.append(Box(-4.0, 4.0, -1.4, 1.4, -100.0, -0.5))        # Bottom_Plate
            box_list.append(Box(-100, 0.5, -1.4, 1.4, -100.0, 100.0))       # Left_Plate_Big
            box_list.append(Box(-100.0, 100.0, 10.4, 100, -100.0, 100.0))   # Back_Plate
            box_list.append(Box(-100.0, 100.0, -100, -10.4, -100.0, 100.0)) # Front_Plate
        else:
            raise ValueError

    elif (denticity == 2) and (input_topology == [2, 2]) and bool_placed_boxes == True: # place bidentate in the plane to the left
        if build_options == 'slab':
            logging.debug("Bidenatate Rotator: Right --> slab chosen")
            box_list.append(Box(-100.0, 100.0, -1.0, 1.0, 0.5, 1.0))        # Top_Plate
            box_list.append(Box(-100.0, 100.0, -1.0, 1.0, -1.0, -0.5))      # Bottom_Plate
            box_list.append(Box(0.1, 100, -1.0, 1.0, -1.0, 1.0))
            box_list.append(Box(0.5, 100, -100, 100, -100, 100))
        elif build_options == 'horseshoe':
            logging.debug("Bidenatate Rotator: Left --> horseshoe chosen")
            box_list.append(Box(-4.0, 4.0, -1.4, 1.4, 0.5, 100))            # Top_Plate
            box_list.append(Box(-4.0, 4.0, -1.4, 1.4, -100.0, -0.5))        # Bottom_Plate
            box_list.append(Box(-0.5, 100, -1.4, 1.4, -100.0, 100.0))       # Right_Plate_Big
            box_list.append(Box(-100.0, 100.0, 10.4, 100, -100.0, 100.0))   # Back_Plate
            box_list.append(Box(-100.0, 100.0, -100, -10.4, -100.0, 100.0)) # Front_Plate
        else:
            raise ValueError


    elif (denticity == 2) and ((input_topology == [3, 2, 0]) or (input_topology == [3, 2, 1])):
        positive_y_edge = 0.1
        negative_y_edge = -0.1
        positive_x_frontier = 0.5
        negative_x_frontier = -3.0
        negative_z_frontier = -3.0
        box_list.append(Box(negative_x_frontier, positive_x_frontier, -100.0, negative_y_edge, negative_z_frontier, 1.0))   # Bottom -y
        box_list.append(Box(negative_x_frontier, positive_x_frontier, positive_y_edge, 100.0, negative_z_frontier, 1.0))    # Bottom +y
        box_list.append(Box(0.0, 2.0, -2.0, 2.0, -1.0, 1.0))                                                                #Top Plate
        """
        box_list.append(Box(-2.0, -1.6, -100.0, -0.5, -5.0, 1.1))   # -x-y
        box_list.append(Box(-2.0, -1.6, 0.5, 100.0, -5.0, 1.1))     # -x+y
        box_list.append(Box(-1.0, 100, -100, 100, -100, 100))       # +x
        """

    return box_list
