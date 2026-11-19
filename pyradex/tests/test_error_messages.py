import pytest

from ..utils import check_unexpected_kwargs


def test_check_unexpected_kwargs_noop():
    check_unexpected_kwargs(lambda a, b: None, {})


def test_check_unexpected_kwargs_alias_and_fuzzy():
    def func(temperature=None, column=None):
        pass

    with pytest.raises(TypeError, match="'tkin' -> 'temperature'"):
        check_unexpected_kwargs(func, {'tkin': 1}, aliases={'tkin': 'temperature'})

    with pytest.raises(TypeError, match="'colum' -> 'column'"):
        check_unexpected_kwargs(func, {'colum': 1})

    with pytest.raises(TypeError, match="Valid keyword arguments are: temperature, column"):
        check_unexpected_kwargs(func, {'zzzzzz': 1})


def _radex():
    pytest.importorskip('pyradex.radex.radex')
    from ..core import Radex
    return Radex


def test_radex_init_unexpected_kwarg():
    Radex = _radex()
    with pytest.raises(TypeError, match="'tkin' -> 'temperature'"):
        Radex(datapath='examples/', species='co', total_density=1e6,
              column=1e16, tkin=100)


def test_radex_init_species_passed_positionally():
    Radex = _radex()
    with pytest.raises(TypeError, match="species='co'"):
        Radex('co', total_density=1e6, column=1e16)


def test_radex_init_conflicting_densities_message():
    Radex = _radex()
    with pytest.raises(ValueError, match="density=.*total_density="):
        Radex(datapath='examples/', species='co', total_density=1e6,
              column=1e16, density=1e3)


def test_radex_set_params_unexpected_kwarg():
    Radex = _radex()
    R = Radex(datapath='examples/', species='co', total_density=1e6,
              column=1e16, temperature=20)
    with pytest.raises(TypeError, match="'tkin' -> 'temperature'"):
        R(tkin=30)
