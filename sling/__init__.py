from pkg_resources import get_distribution


try:
    __version__ = get_distribution('sling').version
except:
    __version__ = 'local'

__all__ = [
    'prepare',
    'scan',
    'filter',
    'group',
    'tasks',
]

from sling import *
