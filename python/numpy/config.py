import os

use_cupy = os.environ.get("USE_CUPY", False)

if use_cupy:
    try:
        import cupy as np
    except (ImportError, ModuleNotFoundError):
        import numpy as np
else:
    import numpy as np
