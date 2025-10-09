from __future__ import annotations


class TransformerChecksMixin:
    def _transform_checks(self, x, /, *fitted_props, **check_params):
        from sklearn.utils.validation import check_is_fitted

        if x is not None:
            x = self._validate_data(x, reset=False, **check_params)
        check_is_fitted(self, *fitted_props)
        return x
