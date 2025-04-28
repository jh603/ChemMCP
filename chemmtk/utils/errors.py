import functools
import inspect

class ChemMTKError(Exception): ...
"""Errors related to the ChemMcpToolkit library."""


class ChemMTKFatalError(ChemMTKError): ...
"""Fatal errors are unrecoverable errors that should not be caught."""

class ChemMTKGeneralError(ChemMTKError): ...
"""General errors are recoverable errors that can be caught."""

class ChemMTKRemoteServerDownError(ChemMTKGeneralError): ...
"""The remote server is down or unreachable."""

class ChemMTKToolInitError(ChemMTKGeneralError): ...

class ChemMTKToolProcessError(ChemMTKGeneralError): ...

class ChemMTKApiNotFoundError(ChemMTKGeneralError): ...

class ChemMTKInputError(ChemMTKGeneralError): ...

class ChemMTKOutputError(ChemMTKGeneralError): ...

class ChemMTKSearchFailError(ChemMTKGeneralError): ...
"""The search on an external resource cannot get reasonable results."""


def catch_errors(func):
    """Wrap sync or async func so that MyError becomes a returned string."""
    if inspect.iscoroutinefunction(func):
        @functools.wraps(func)
        async def wrapper(*args, **kwargs):
            try:
                return await func(*args, **kwargs)
            except ChemMTKError as e:
                return f"Error ({e.__class__.__name__}): {e}"
        return wrapper
    else:
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            try:
                return func(*args, **kwargs)
            except ChemMTKError as e:
                return f"Error ({e.__class__.__name__}): {e}"
        return wrapper
