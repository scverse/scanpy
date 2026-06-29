{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. add toctree option to make autodoc generate the pages

..  autoclass:: {{ objname }}

    {% block attributes %}
    {%- for item in attributes %}
    {%- if loop.first %}
    .. rubric:: Attributes

    ..  autosummary::
        :toctree: .
    {% endif %}
        ~{{ name }}.{{ item }}
    {%- endfor %}
    {% endblock %}

    {% block methods %}
    {%- for item in methods if item != '__init__' %}
    {%- if loop.first %}
    .. rubric:: Methods

    ..  autosummary::
        :toctree: .
    {% endif %}
        ~{{ name }}.{{ item }}
    {%- endfor %}
    {% endblock %}
