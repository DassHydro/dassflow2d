{{ fullname }}
{{ underline }}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}

{% for item in members %}
{%- if item in ['msh'] %}
	.. autoclass :: {{ module }}.{{ objname }}.{{ item }}

{%- endif -%}
{%- endfor %}

..
   HACK -- meth the point here is that we don't want this to appear in the output, but the autosummary should still generate the pages.
   .. autosummary::
      :toctree:
      :noindex:
      {% for item in members %}
      {%- if not item.startswith('_') %}
      {%- if not item in ['msh'] %}
      {{ name }}.{{ item }}
      {%- endif -%}
      {%- endif -%}
      {%- endfor %}

	
      
      
