from types import MappingProxyType
from typing import Any, Mapping, Sequence, NamedTuple

from docutils import nodes
from docutils.parsers.rst.directives import class_option
from docutils.parsers.rst.states import Inliner
from sphinx.application import Sphinx
from sphinx.config import Config


class AutoLink(NamedTuple):
    url_template: str
    title_template: str = '{}'
    options: Mapping[str, Any] = MappingProxyType({'class': class_option})

    def __call__(
        self,
        name: str,
        rawtext: str,
        text: str,
        lineno: int,
        inliner: Inliner,
        options: Mapping[str, Any] = MappingProxyType({}),
        content: Sequence[str] = (),
    ):
        url = self.url_template.format(text)
        title = self.title_template.format(text)
        node = nodes.reference(rawtext, title, refuri=url, **options)
        return [node], []


def register_links(app: Sphinx, config: Config):
    gh_url = 'https://github.com/{github_user}/{github_repo}'.format_map(
        config.html_context
    )
    app.add_role('pr', AutoLink(f'{gh_url}/pull/{{}}', 'PR {}'))
    app.add_role('issue', AutoLink(f'{gh_url}/issues/{{}}', 'issue {}'))


def setup(app: Sphinx):
    app.connect('config-inited', register_links)
