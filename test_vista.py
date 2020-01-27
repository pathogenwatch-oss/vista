from unittest import TestCase

from vista.vista_utils import overlaps


class Test(TestCase):
    def test_overlaps(self):
        self.assertTrue(overlaps((1, 100), (20, 120), 50))
        self.assertFalse(overlaps((1, 100), (20, 120), 100))
        self.assertFalse(overlaps((1, 110), (120, 200), 60))
        self.assertTrue(overlaps((1, 100), (20, 80), 40))
