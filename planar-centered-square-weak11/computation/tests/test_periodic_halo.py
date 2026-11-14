import importlib.util
from fractions import Fraction as F
from pathlib import Path
import unittest

spec = importlib.util.spec_from_file_location('periodic_halo', Path(__file__).parents[1] / 'periodic_halo.py')
h = importlib.util.module_from_spec(spec)
spec.loader.exec_module(h)

class HaloTests(unittest.TestCase):
    def test_union_overlap_touch_and_degenerate_input(self):
        self.assertEqual(h.union_area([]), 0)
        self.assertEqual(h.union_area([(0,2,0,2),(1,3,1,3),(3,4,0,1)]), 8)

    def test_one_atom_period_border_and_square_root_equality(self):
        data={'Px':'10','Py':'10','lam':'1','pts':[['0','0','1']]}
        r=h.calculate(data,bits=12,tail=F(5))
        self.assertLess(F(r['inner_area']), 1)
        self.assertGreater(F(r['outer_area']), 1)
        self.assertLess(F(r['outer_area'])-F(r['inner_area']), F(1,100))
        self.assertGreater(r['counts']['root_checked'],0)
        shifted={'Px':'10','Py':'10','lam':'1','pts':[['10','-10','1']]}
        r2=h.calculate(shifted,bits=12,tail=F(5))
        self.assertEqual((r['inner_area'],r['outer_area']), (r2['inner_area'],r2['outer_area']))

    def test_collision_and_symmetry(self):
        base={'Px':'10','Py':'10','lam':'1','pts':[['0','0','1'],['1','2','1/2']]}
        collision={'Px':'10','Py':'10','lam':'1','pts':[['10','0','1/2'],['0','0','1/2'],['1','2','1/2']]}
        flipped={'Px':'10','Py':'10','lam':'1','pts':[['0','0','1'],['2','1','1/2']]}
        results=[h.calculate(d,12,F(5)) for d in (base,collision,flipped)]
        self.assertEqual(len({(r['inner_area'],r['outer_area']) for r in results}),1)

    def test_zero_mass_and_invalid_inputs(self):
        base={'Px':'10','Py':'10','lam':'1','pts':[['0','0','1']]}
        zero={'Px':'10','Py':'10','lam':'1','pts':[['0','0','1'],['3','4','0']]}
        a,b=(h.calculate(d,12,F(5)) for d in (base,zero))
        self.assertEqual((a['inner_area'],a['outer_area']), (b['inner_area'],b['outer_area']))
        for d in ({'Px':'0','Py':'10','lam':'1','pts':[['0','0','1']]},
                  {'Px':'10','Py':'10','lam':'1','pts':[['0','0','-1']]},
                  {'Px':'10','Py':'10','lam':'1','pts':[['0','0','0']]}):
            with self.assertRaises(ValueError): h.calculate(d,12,F(5))

    def test_scaling_and_nonunit_level(self):
        base={'Px':'10','Py':'10','lam':'1','pts':[['0','0','1'],['1','2','1/2']]}
        scaled={'Px':'20','Py':'20','lam':'1','pts':[['0','0','4'],['2','4','2']]}
        level={'Px':'10','Py':'10','lam':'3','pts':[['0','0','3'],['1','2','3/2']]}
        a=h.calculate(base,12,F(5))
        b=h.calculate(scaled,11,F(10))
        c=h.calculate(level,12,F(5))
        for key in ('inner_area','outer_area'):
            self.assertEqual(F(b[key]),4*F(a[key]))
            self.assertEqual(F(c[key]),F(a[key]))
        for key in ('inner_ratio','outer_ratio'):
            self.assertEqual(F(a[key]),F(b[key]))
            self.assertEqual(F(a[key]),F(c[key]))

    def test_tail_must_be_certified(self):
        with self.assertRaises(ValueError):
            h.calculate({'Px':'1','Py':'1','lam':'1','pts':[['0','0','1']]},12,F(5))

    def test_cap_and_empty_inner(self):
        d={'Px':'10','Py':'10','lam':'1','pts':[['5','5','1/1048576']]}
        r=h.calculate(d,2,F(1))
        self.assertEqual(F(r['inner_area']), 0)
        self.assertGreater(F(r['outer_area']), 0)
        self.assertGreater(r['counts']['inner_empty'],0)
        self.assertEqual(r['counts']['outer_empty'],0)

if __name__ == '__main__': unittest.main()
