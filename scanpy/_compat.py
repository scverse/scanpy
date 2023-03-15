from packaging import version


def pkg_metadata(package):
    from importlib.metadata import metadata as m

    return m(package)


def pkg_version(package):
    from importlib.metadata import version as v

    return version.parse(v(package))
