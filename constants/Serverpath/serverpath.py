"""
In progress, hopefully this is going to become a serverpath anytime soon
"""
from pathlib import Path


def get_serverpath(possible_serverpaths: list):
    """
    # Returns serverpath on this machine based on trying a list of paths and seeing which path exists.
    """
    serverpath = None

    for path in possible_serverpaths:
        if Path(path).exists():
            serverpath = path
            break

    if serverpath is None:
        raise FileNotFoundError(f'No existing serverpath found in {possible_serverpaths}. Please doublecheck and add the path to your local copy of the ccem server to the list of possible serverpaths.')

    return serverpath


possible_serverpaths = [
                        "/Users/felixk/Documents/PhD/Data/CreateTMC",
                        '/Users/timosommer/PhD/projects/RCA/projects/CreateTMC/ccem_server'
                        ]
serverpath = get_serverpath(possible_serverpaths=possible_serverpaths)