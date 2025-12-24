{{ fullname | escape | underline}}

.. currentmodule:: {{ module }}

.. add toctree option to make autodoc generate the pages

.. autoclass:: {{ objname }}

{% block attributes %}
{% if attributes %}
Attributes table
~~~~~~~~~~~~~~~~

.. autosummary::
{% for item in attributes %}
    ~{{ name }}.{{ item }}
{%- endfor %}
{% endif %}
{% endblock %}

{% block methods %}
{% if methods %}
Methods table
~~~~~~~~~~~~~

.. autosummary::
{% for item in methods %}
    {%- if item != '__init__' %}
    ~{{ name }}.{{ item }}
    {%- endif -%}
{%- endfor %}
{% endif %}
{% endblock %}

{% block attributes_documentation %}
{% if attributes %}
Attributes
~~~~~~~~~~

{% for item in attributes %}

.. autoattribute:: {{ [objname, item] | join(".") }}
{%- endfor %}

{% endif %}
{% endblock %}

{% block methods_documentation %}
{% if methods %}
Methods
~~~~~~~

{% for item in methods %}
{%- if item != '__init__' %}

.. automethod:: {{ [objname, item] | join(".") }}
{%- endif -%}
{%- endfor %}

{% endif %}
{% endblock %}
