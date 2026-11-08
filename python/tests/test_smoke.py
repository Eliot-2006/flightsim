def test_scientific_stack_and_jsbsim_import():
    import jsbsim
    import numpy as np

    assert np.dot([1.0, 0.0, 0.0], [1.0, 0.0, 0.0]) == 1.0
    assert jsbsim.FGFDMExec is not None
