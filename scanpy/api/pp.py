# we need this dummy module as otherwise, sphinx cannot produce the overview of
# the api module
#
# importing all ("*") is ok, as scanpy.preprocessing.__init__ is carefully maintained
# to only contain the functions of the api
from ..preprocessing import *
