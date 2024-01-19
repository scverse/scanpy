from __future__ import annotations


class TransformerChecksMixin:
    def _transform_checks(self, X, *fitted_props, **check_params):
        from sklearn.utils.validation import check_is_fitted

        if X is not None:
            X = self._validate_data(X, reset=False, **check_params)
        check_is_fitted(self, *fitted_props)
        return X
