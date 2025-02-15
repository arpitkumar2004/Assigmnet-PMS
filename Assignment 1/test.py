import unittest
import numpy as np
from C_b import Distillation

class TestDistillation(unittest.TestCase):

    def setUp(self):
        self.distillation = Distillation()

    def test_calculate_Pij(self):
        C1 = np.array([81.768, 74.475, 88.134])
        C2 = np.array([-6876, -7164.3, -8438.6])
        C3 = np.array([-8.7078, -7.327, -9.0766])
        C4 = np.array([7.1926e-6, 3.1340e-6, 8.3303e-18])
        C5 = np.array([2, 2, 6])
        T_j = np.array([65, 98, 97])
        Pij_sat = self.distillation.calculate_Pij(C1, C2, C3, C4, C5, T_j)
        self.assertEqual(Pij_sat.shape, (3, 3))

    def test_calculatate_Kij(self):
        Pij_sat = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        Kij = self.distillation.calculatate_Kij(Pij_sat)
        self.assertEqual(Kij.shape, (3, 3))

    def test_calculate_Lj(self):
        feed_tray = 2
        Vj = np.array([0, 500, 500])
        Lj = self.distillation.calculate_Lj(feed_tray, Vj)
        self.assertEqual(Lj.shape, (3,))

    def test_calculate_Aj(self):
        Lj = np.array([100, 200, 300])
        Aj = self.distillation.calculate_Aj(Lj)
        self.assertEqual(Aj.shape, (3,))

    def test_calculate_Bij(self):
        Lj = np.array([100, 200, 300])
        Kij = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        Bij = self.distillation.calculate_Bij(Lj, Kij)
        self.assertEqual(Bij.shape, (3, 3))

    def test_calculate_Cij(self):
        Vj = np.array([100, 200, 300])
        Kij = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        Cij = self.distillation.calculate_Cij(Vj, Kij)
        self.assertEqual(Cij.shape, (3, 3))

    def test_calculate_Dij(self):
        feed_tray = 2
        Dij = self.distillation.calculate_Dij(feed_tray)
        self.assertEqual(Dij.shape, (3, 3))

    def test_thomas_algorithm(self):
        A = np.array([0, 1, 1])
        B = np.array([4, 4, 4])
        C = np.array([1, 1, 0])
        D = np.array([7, 8, 9])
        x = self.distillation.thomas_algorithm(A, B, C, D)
        self.assertEqual(x.shape, (3,))

    def test_check_summation(self):
        Xij = np.array([0.6, 0.2, 0.2])
        result = self.distillation.check_summation(Xij)
        self.assertTrue(result)

    def test_calculate_Tj(self):
        A = np.array([5.20409, 5.3229, 5.31384])
        B = np.array([1581.341, 1670.409, 1690.864])
        C = np.array([-33.5, -40.191, -51.804])
        Xij = np.array([[0.6, 0.2, 0.2], [0.5, 0.3, 0.2], [0.4, 0.4, 0.2]])
        Tj = self.distillation.calculate_Tj(A, B, C, Xij)
        self.assertEqual(Tj.shape, (3,))

    def test_calculate_Vj_heat_balance(self):
        Lj = np.array([100, 200, 300])
        hL = np.array([10, 20, 30])
        hV = np.array([40, 50, 60])
        hF = 70
        Q = 0
        feed_tray = 2
        Vj = self.distillation.calculate_Vj_heat_balance(Lj, hL, hV, hF, Q, feed_tray)
        self.assertEqual(Vj.shape, (3,))

    def test_calculate_enthalpies(self):
        Xij = np.array([[0.6, 0.2, 0.2], [0.5, 0.3, 0.2], [0.4, 0.4, 0.2]])
        Kij = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        hij = np.array([-238.8, -277.6, -326.6])
        hL, hV, hF = self.distillation.calculate_enthalpies(Xij, Kij, hij)
        self.assertEqual(hL.shape, (3,))
        self.assertEqual(hV.shape, (3,))
        self.assertIsInstance(hF, float)

    def test_check_convergence(self):
        Vj = np.array([100, 200, 300])
        Tj = np.array([65, 98, 97])
        Vj_prev = np.array([100, 200, 300])
        Tj_prev = np.array([65, 98, 97])
        result = self.distillation.check_convergence(Vj, Tj, Vj_prev, Tj_prev)
        self.assertTrue(result)

if __name__ == '__main__':
    unittest.main()