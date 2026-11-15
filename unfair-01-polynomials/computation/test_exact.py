import unittest
from algebraic_verify import RealField, verify
from support_search import classify
from verify_support_ledger import verify_record


class ExactChecks(unittest.TestCase):
    def certificate(self,A,B):
        return {'minpoly':[0,1],'interval':[-1,1],'A':[[x] for x in A],'B':[[x] for x in B]}

    def test_valid_and_probability_return(self):
        r=verify(self.certificate([1,1],[1,0,1]))
        self.assertFalse(r['counterexample'])
        self.assertEqual(r['uniform_sum_probability'],'1/4')
        self.assertEqual(r['sum_support'],[0,1,2,3])

    def test_invalid_products(self):
        for A,B in [([1,1],[1,1]),([1,-1,1],[1,1]),([2,2],['1/2',0,'1/2']),([1,'999999999/1000000000',1],[1,0,0,1])]:
            with self.assertRaises(ValueError):
                verify(self.certificate(A,B))

    def test_no_float(self):
        with self.assertRaises(ValueError):
            verify(self.certificate([1,0.5,1],[1,1]))

    def test_algebraic_sign_and_identity(self):
        f=RealField([-2,0,1],[1,2])
        t=f.norm([0,1])
        self.assertEqual(f.mul(t,t),f.norm([2]))
        self.assertEqual(f.sign(f.norm([-1,1])),1)
        self.assertEqual(f.sign(f.norm([1,-1])),-1)
        self.assertEqual(f.sign(f.norm([-2,0,1])),0)
        self.assertEqual(f.mul(t,f.inverse(t)),f.one)
        with self.assertRaises(ValueError):
            RealField([-1,0,1],[0,2])
        with self.assertRaises(ValueError):
            RealField([-2,0,1],[-2,2])

    def test_algebraic_false_candidate(self):
        with self.assertRaises(ValueError):
            verify({'minpoly':[-2,0,1],'interval':[1,2],
                    'A':[[1],[-1,1],[1]],'B':[[1],[2,-1],[1]]})

    def test_singleton_certificates(self):
        for I,J in [((0,1),(0,2)),((0,1),(0,1)),((0,2,5),(0,2,3,5,7))]:
            status,witness=classify(I,J)
            verify_record([I,J,status,witness])
        with self.assertRaises(AssertionError):
            verify_record([[0,1],[0,2],'excluded',[1,1,0]])


if __name__=='__main__':
    unittest.main()
