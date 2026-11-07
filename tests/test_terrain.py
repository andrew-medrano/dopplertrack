"""Terrain tiles are decoded once: the PNG becomes a float32 .npy that later reads memory-map."""
import numpy as np
from PIL import Image

from dopplertrack import terrain


def _fake_download(self, x, y):
    p = self._path(x, y)
    if p.exists():
        return
    elev = 120.0 + 0.5 * np.arange(256)[None, :] + np.zeros((256, 1))       # rises 0.5 m per pixel eastwards
    v = elev + 32768
    rgb = np.stack([np.floor(v / 256), np.floor(v % 256), np.round((v % 1) * 256)], -1).astype(np.uint8)
    Image.fromarray(rgb, "RGB").save(p)


def test_tiles_are_decoded_once(tmp_path, monkeypatch):
    monkeypatch.setattr(terrain.Terrain, "_download", _fake_download)
    t = terrain.Terrain(zoom=9, cache_dir=tmp_path)
    lat, lon = np.array([40.0, 40.0]), np.array([-73.6, -73.5])   # one tile at zoom 9
    e = t.elevation(lat, lon)
    assert np.all((e >= 120.0) & (e <= 248.0)) and e[1] > e[0]
    npys = list(t.cache.glob("*.npy"))
    assert len(npys) == 1 and np.load(npys[0]).dtype == np.float32
    for p in t.cache.glob("*.png"):
        p.unlink()
    monkeypatch.setattr(terrain.Terrain, "_download", lambda self, x, y: (_ for _ in ()).throw(AssertionError("downloaded")))
    again = terrain.Terrain(zoom=9, cache_dir=tmp_path)
    assert np.allclose(again.elevation(lat, lon), e)
    assert isinstance(again._mem[next(iter(again._mem))], np.memmap)
