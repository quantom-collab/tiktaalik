import numpy as np
from tiktaalik.model import Hu, Hd, Hs, Hg

# GOLDEN MASTER: frozen snapshot of the model output at fixed kinematics.
# Values captured from the built tiktaalik (master). This is the one absolute-value
# anchor the valence+sea / singlet relations can't give (those are relative and would
# miss a global rescale).
X = np.array([0.15, 0.3, 0.5, 0.7, 0.85])
XI = 0.1
T = 10.0

_GOLDEN = {
    'Hu': np.array([10109969404.669182, 211779.11790729326, 387.05394531120265, 3.3156034556609213, 0.0671045335440013]),
    'Hd': np.array([5779124133.685159, 92975.34213670118, 110.09292521835022, 0.4999832650387941, 0.003948734991152491]),
    'Hs': np.array([7.221888263975643, 0.31442024228090304, 0.013379863055888563, 0.00031769368797540665, 3.493041731417179e-06]),
    'Hg': np.array([20.23639816694913, 2.0717825860192534, 0.3335970756597817, 0.02985204989090562, 0.0010296208783984184]),
}


def test_model_golden_master():
    for name, fn in [('Hu', Hu), ('Hd', Hd), ('Hs', Hs), ('Hg', Hg)]:
        got = np.ravel(np.asarray(fn(X, XI, T)))
        assert np.allclose(got, _GOLDEN[name], rtol=1e-12, atol=1e-14), name
