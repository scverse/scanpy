__all__ = ['urlretrieve']

try:
    from urllib.request import urlretrieve
except ImportError:  # Python 2
    from urllib import urlretrieve
