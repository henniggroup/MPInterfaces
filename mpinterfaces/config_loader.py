import os
from monty.serialization import loadfn

def load_configurations():
    expected_keys = [
        'PMG_MAPI_KEY', 'PMG_VASP_PSP_DIR', 'username', 'normal_binary',
        'twod_binary', 'vdw_kernel', 'potentials', 'queue_system', 'queue_template'
    ]
    configurations = {}

    # Attempt to load from .pmgrc.yaml
    pmgrc_path = os.path.join(os.path.expanduser('~'), '.pmgrc.yaml')
    if os.path.isfile(pmgrc_path):
        try:
            pmgrc_config = loadfn(pmgrc_path)
            for key in expected_keys:
                configurations[key] = pmgrc_config.get(key)
        except Exception as e:
            print(f"Error loading configurations from {pmgrc_path}: {e}")

    # Override with environment variables or set to None if not found
    for key in expected_keys:
        configurations[key] = os.getenv(key.upper(), configurations.get(key))

    # Ensure all expected keys are present, even if as None
    for key in expected_keys:
        if key not in configurations:
            configurations[key] = None

    return configurations

CONFIG = load_configurations()
