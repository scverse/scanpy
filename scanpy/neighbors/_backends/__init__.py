from types import MappingProxyType

from . import third_party

BACKENDS = MappingProxyType(third_party.BACKENDS)
ALGORITHMS = MappingProxyType(third_party.BACKENDS)
