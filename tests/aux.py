from pathlib import Path

def getpath(*a):
    # Package root
    d = Path(__file__).parent.resolve()
    return d.joinpath(*a)

