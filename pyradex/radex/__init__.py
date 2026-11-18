"""pyradex.radex package initializer.

Try to import the compiled extension module `pyradex.radex.radex` (the
extension shared object) in a way that avoids recursively re-importing this
package's __init__ module. If the compiled extension is missing or fails to
import, raise a clear ImportError so the user can rebuild the extension.
"""

try:
    from . import radex
except Exception as _exc:  # pragma: no cover - runtime import failures
	raise ImportError(
		"Could not import the compiled radex extension (pyradex.radex.radex). "
		"Build the extension with `python setup.py install_myradex` or see "
		"install_radex.py for instructions. Original error: %r" % _exc
	)
