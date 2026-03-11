import os
import shutil
import pytest
import pandas as pd
import numpy as np
import phydrus as ps

@pytest.fixture
def test_workspace(tmpdir):
    ws = str(tmpdir)
    yield ws
    shutil.rmtree(ws, ignore_errors=True)

def test_write_input_no_nans(test_workspace):
    """
    Test that write_input() handles NaNs correctly and produces files
    without 'NaN' strings which break Fortran reading.
    """
    # Dynamically find the hydrus binary in the project root
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    hydrus_exe = os.path.join(project_root, "hydrus1d", "bin", "hydrus")
    
    # Create simple model
    ml = ps.Model(exe_name=hydrus_exe, ws_name=test_workspace,
                  name="hydrus", description="Test Model")
    
    # Add basic atmospheric info with some NaNs
    atmosphere = pd.DataFrame({
        "tAtm": [1.0, 2.0],
        "Prec": [0.0, 1.0],
        "rSoil": [np.nan, 0.05],
        "rRoot": [pd.NA, np.nan],
        "hCritA": [100000, 100000],
        "rB": [0, 0],
        "hB": [0, 0],
        "ht": [0, 0],
        "tTop": [0, 0],
        "tBot": [0, 0],
        "Ampl": [0, 0],
        "cTop": [0.0, 0.0],
        "cBot": [0.0, pd.NA]
    })
    ml.add_atmospheric_bc(atmosphere)
    ml.add_waterflow(model=0, hysteresis=0)
    
    mat_df = pd.DataFrame({"thr": [0.0], "ths": [0.338], "Alfa": [0.011], 
                           "n": [1.47], "Ks": [13.0], "l": [0.5]}, index=[1])
    ml.add_material(mat_df)
    ml.add_time_info(tinit=0.0, tmax=2.0, dt=0.01)
    
    # Write input files
    ml.write_selector()
    ml.write_atmosphere()
    
    # Assert ATMOSPH.IN does not contain 'NaN' strings
    atm_path = os.path.join(test_workspace, "ATMOSPH.IN")
    assert os.path.exists(atm_path)
    
    with open(atm_path, "r") as f:
        content = f.read()
        
    # Pandas NaNs must be explicitly removed or cast to 0
    assert "NaN" not in content, "Found 'NaN' string in ATMOSPH.IN! Fortran will fail to read this."
    assert "<NA>" not in content, "Found '<NA>' string in ATMOSPH.IN! Fortran will fail to read this."
    
    # Check that EOF has proper formatting and newline
    assert "end*** END OF INPUT FILE ATMOSPH.IN" in content

def test_write_profile_no_nans(test_workspace):
    project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
    hydrus_exe = os.path.join(project_root, "hydrus1d", "bin", "hydrus")
    ml = ps.Model(exe_name=hydrus_exe, ws_name=test_workspace, name="hydrus")
    ml.add_waterflow(model=0, hysteresis=0)
    
    profile = pd.DataFrame({
        "x": [0, -10, -20],
        "h": [-100, np.nan, -100],
        "Mat": [1, 1, np.nan],
        "Lay": [1, 1, 1],
        "CosAlf": [1.0, 1.0, 1.0],
        "W": [pd.NA, 1.0, 1.0]
    })
    
    ml.add_profile(profile)
    ml.write_profile()
    
    prof_path = os.path.join(test_workspace, "PROFILE.DAT")
    assert os.path.exists(prof_path)
    
    with open(prof_path, "r") as f:
        content = f.read()
        
    assert "NaN" not in content, "Found 'NaN' in PROFILE.DAT!"
    assert "<NA>" not in content, "Found '<NA>' in PROFILE.DAT!"
