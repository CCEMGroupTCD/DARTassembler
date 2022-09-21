
def coordinates_to_xyz_str(coordinates: dict):
    """
    returns a string that can be written into an .xyz file
    """
    str_ = f"{len(list(coordinates.keys()))}\n \n"
    for coord in coordinates.values():
        str_ += f"{coord[0]} \t {coord[1][0]} \t {coord[1][1]} \t {coord[1][2]} \n"

    return str_


