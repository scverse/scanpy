from sklearn.utils import deprecated

from ._mnn_correct import mnn_correct
from ._bbknn import bbknn
from ._dca import dca
from ._harmony_integrate import harmony_integrate
from ._magic import magic
from ._scanorama_integrate import scanorama_integrate
from ._hashsolo import hashsolo
from ...preprocessing import _scrublet

scrublet = deprecated("Import from sc.pp instead")(_scrublet.scrublet)
scrublet_simulate_doublets = deprecated("Import from sc.pp instead")(
    _scrublet.scrublet_simulate_doublets
)
