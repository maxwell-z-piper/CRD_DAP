"""Recognize verified presentation-only Script 07 upgrades."""

from .bl_population_io import fingerprint

# Pin the exact driver/exporter pairs that left the fits and draw checkpoints
# compatible. All other code and scientific settings must still match.
_COMPATIBLE_REVISIONS = (
    ('961670255f6f88804a55c337fe1e921a82883a03b487c0e87fd8db3bf33f9ed1',
     'acaac550ba190a5a022d46723d9f83998f8b604a172fe8251b882063bfb86a55'),
    ('0ef5be2a2e91bbeea59cca9a36f6d8e283c345ade7a68174ad9ebff3afdc036d',
     'a50411117f919bfb1f945f0068dd80e102868f2448d3536b6307e53a5ac92acd'),
    ('0ef5be2a2e91bbeea59cca9a36f6d8e283c345ade7a68174ad9ebff3afdc036d',
     '5681ee86b238c310e0423945986d8effa622710b55d66e1ccbb559d1182d5e01'),
    ('0ef5be2a2e91bbeea59cca9a36f6d8e283c345ade7a68174ad9ebff3afdc036d',
     'ebbdd2e8cb7efffcca9886a10e155bd8128b7c7225450cad55265be911141bd6'),
)
_PLOT_REVISION = (
    '0ef5be2a2e91bbeea59cca9a36f6d8e283c345ade7a68174ad9ebff3afdc036d',
    'b0491c1a368393854be995100c224564e14334ce41e788db8f63888b4797518a',
)


def matches_checkpoint(saved, current):
    if saved['fingerprint'] == current['fingerprint']:
        return True
    contract = current['contract']
    code = contract['code']
    if (contract['driver_sha256'], code['script07_products.py']) != _PLOT_REVISION:
        return False
    for driver, exporter in _COMPATIBLE_REVISIONS:
        previous_contract = dict(contract, driver_sha256=driver,
                                 code=code | {'script07_products.py': exporter})
        if saved['fingerprint'] == fingerprint(previous_contract):
            return True
    return False
