{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

..  autoclass:: {{ objname }}
    :show-inheritance:

    {% block attributes %}
    .. rubric:: Members
    {% for item in attributes if item[0] is upper %}
    .. autoattribute:: {{ item }}
    {%- endfor %}

    {%- for item in attributes if item[0] is not upper %}
    {%- if loop.first %}
    .. rubric:: Attributes
    ..  autosummary::
        :toctree: .
    {% endif %}
        ~{{ name }}.{{ item }}
    {%- endfor %}
    {% endblock %}

    {% block methods %}
    {%- for item in methods if item != "__init__" and item not in inherited_members %}
    {%- if loop.first %}
    .. rubric:: Methods

    ..  autosummary::
        :toctree: .
    {% endif %}
        ~{{ name }}.{{ item }}
    {%- endfor %}
    {% endblock %}
