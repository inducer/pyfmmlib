from importlib import metadata


VERSION_TEXT = metadata.version("pyfmmlib")
VERSION = tuple([int(v) for v in VERSION_TEXT.split(".")])
VERSION_STATUS = ""
