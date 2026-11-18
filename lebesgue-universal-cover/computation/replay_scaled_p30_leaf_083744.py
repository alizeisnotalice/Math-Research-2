"""Check closed prefixes of one original 0.83744 leaf with contracted material."""
from hashlib import sha256
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import replay_p30_prefix_180_24225 as base
import tomar_p30_mesh_guard as mesh
from tomar_p30_scaled_material import scaled_material

if __name__ == '__main__':
    for name, digest in (
        ('replay_p30_prefix_180_24225.py', '50c1340be3d58c1ae4f2f213bb4315915d8935ad379742221064222141389941'),
        ('tomar_p30_mesh_guard.py', '6781862477e5fb7fd8cadc63e91a1c34c0a530508bf3fe785e95efc90a4c95d7'),
        ('tomar_p30_scaled_material.py', 'da3fdc4472aaf0b09e020b6a3a7b467f227fe695c84ad0c4f502244e0f20cc2b'),
    ):
        assert sha256((HERE / name).read_bytes()).hexdigest() == digest
    selection = json.loads(Path(sys.argv[-1]).read_bytes())
    assert selection.get('terminal_guard') == 'mesh64-dense'
    base.PART, base.POS = selection['part'], selection['pos']
    mesh.dense_material = scaled_material
    base.main()
