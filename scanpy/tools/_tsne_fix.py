# Author: David DeTomaso (https://github.com/deto)
"""\
Fix for sklearn.tsne gradient descent.

This module fixes it by patching the original function with this
modified version. Patch is only applied for versions earlier than 0.19.

Code courtesy of David DeTomaso; available from
https://github.com/YosefLab/FastProject/blob/stable/FastProject/_tsne_fix.py
"""
from types import MappingProxyType
from typing import Callable, Tuple, Optional, Mapping, Any, Iterable

import numpy as np
from scipy import linalg
import sklearn


def _gradient_descent(
    objective: Callable[..., Tuple[int, np.ndarray]],
    p0: np.ndarray,
    it: int,
    n_iter: int,
    objective_error: Optional[Callable[..., float]] = None,
    n_iter_check: int = 1,
    n_iter_without_progress: int = 50,
    momentum: float = 0.5,
    learning_rate: float = 1000.0,
    min_gain: float = 0.01,
    min_grad_norm: float = 1e-7,
    min_error_diff: float = 1e-7,
    verbose: int = 0,
    args: Iterable[Any] = (),
    kwargs: Mapping[str, Any] = MappingProxyType({}),
) -> Tuple[np.ndarray, float, int]:
    """\
    Batch gradient descent with momentum and individual gains.

    Parameters
    ----------
    objective
        Should return a tuple of cost and gradient for a given parameter
        vector. When expensive to compute, the cost can optionally
        be None and can be computed every n_iter_check steps using
        the objective_error function.
    p0
        Initial parameter vector. shape (n_params,)
    it
        Current number of iterations (this function will be called more than
        once during the optimization).
    n_iter
        Maximum number of gradient descent iterations.
    objective_error
        Should return error for a given parameter vector.
    n_iter_check
        Number of iterations before evaluating the global error. If the error
        is sufficiently low, we abort the optimization.
    n_iter_without_progress
        Maximum number of iterations without progress before we abort the
        optimization.
    momentum
        The momentum generates a weight for previous gradients that decays
        exponentially. within (0.0, 1.0)
    learning_rate
        The learning rate should be extremely high for t-SNE! Values in the
        range [100.0, 1000.0] are common.
    min_gain
        Minimum individual gain for each parameter.
    min_grad_norm
        If the gradient norm is below this threshold, the optimization will
        be aborted.
    min_error_diff
        If the absolute difference of two successive cost function values
        is below this threshold, the optimization will be aborted.
    verbose
        Verbosity level.
    args
        Arguments to pass to objective function.
    kwargs
        Keyword arguments to pass to objective function.
    Returns
    -------
    p
        Optimum parameters. shape (n_params,)
    error
        Optimum.
    i
        Last iteration.
    """

    p = p0.copy().ravel()
    update = np.zeros_like(p)
    gains = np.ones_like(p)
    error = np.finfo(np.float).max
    best_error = np.finfo(np.float).max
    best_iter = 0

    for i in range(it, n_iter):
        new_error, grad = objective(p, *args, **kwargs)
        grad_norm = linalg.norm(grad)

        inc = update * grad < 0.0
        dec = np.invert(inc)
        gains[inc] += 0.2
        gains[dec] *= 0.8
        np.clip(gains, min_gain, np.inf, out=gains)
        grad *= gains
        update = momentum * update - learning_rate * grad
        p += update

        if (i + 1) % n_iter_check == 0:
            if new_error is None:
                new_error = objective_error(p, *args)
            error_diff = np.abs(new_error - error)
            error = new_error

            if verbose >= 2:
                m = "[t-SNE] Iteration %d: error = %.7f, gradient norm = %.7f"
                print(m % (i + 1, error, grad_norm))

            if error < best_error:
                best_error = error
                best_iter = i
            elif i - best_iter > n_iter_without_progress:
                if verbose >= 2:
                    print("[t-SNE] Iteration %d: did not make any progress "
                          "during the last %d episodes. Finished."
                          % (i + 1, n_iter_without_progress))
                break
            if grad_norm <= min_grad_norm:
                if verbose >= 2:
                    print("[t-SNE] Iteration %d: gradient norm %f. Finished."
                          % (i + 1, grad_norm))
                break
            if error_diff <= min_error_diff:
                if verbose >= 2:
                    m = "[t-SNE] Iteration %d: error difference %f. Finished."
                    print(m % (i + 1, error_diff))
                break

        if new_error is not None:
            error = new_error

    return p, error, i


sk_ver = []
for c in sklearn.__version__.split("."):
    try:
        ic = int(c)
        sk_ver.append(ic)
    except ValueError:
        pass
sk_ver = tuple(sk_ver)

if sk_ver < (0, 19, 0):
    from sklearn.manifold import t_sne
    t_sne._gradient_descent = _gradient_descent
