import struct
import unittest

import numpy as np

from iplan_dicom.images import ImageVolume
from iplan_dicom.structures import decode_voxel_mask, mask_to_contours
from iplan_dicom.xbrain import FormatError


def encoded(planes):
    data = b'\x04\x08' + struct.pack('<I', 16) + np.eye(4, dtype='<f8').tobytes() + struct.pack('<I', len(planes))
    for runs in planes:
        data += struct.pack('<I', len(runs)) + np.array(runs, dtype='<u4').tobytes()
    return data


class MaskTests(unittest.TestCase):
    def test_alternating_runs_and_empty_plane(self):
        mask = decode_voxel_mask(encoded([[1, 2, 3], [6]]), base_width=3, base_height=2, xy_scale=1)
        np.testing.assert_array_equal(mask.plane(0), [[False, True, True], [False, False, False]])
        self.assertFalse(mask.plane(1).any())
        self.assertEqual(mask.voxel_count, 2)

    def test_invalid_sum_truncation_and_tail_fail(self):
        for data in (encoded([[5]]), encoded([[6]])[:-1], encoded([[6]]) + b'\0'):
            with self.assertRaises(FormatError):
                decode_voxel_mask(data, base_width=3, base_height=2, xy_scale=1)

    def test_scale_and_physical_contour_position(self):
        # 4x4 finer mask, one voxel in row1,col2, native base image2x2.
        mask = decode_voxel_mask(encoded([[6, 1, 9]]), base_width=2, base_height=2, xy_scale=2)
        affine = np.diag([-2., 4., -3., 1.])
        affine[:3, 3] = [10, 20, 30]
        image = ImageVolume(np.zeros((1, 2, 2), dtype='u2'), affine, np.eye(4), (), 'CT', 1, -1024, '')
        contours = mask_to_contours(mask, image)
        self.assertEqual(len(contours), 1)
        z, points = contours[0]
        self.assertEqual(z, 0)
        np.testing.assert_allclose(points.mean(axis=0), [8.5, 21, 30])
        np.testing.assert_allclose(points[:, 0].max() - points[:, 0].min(), 1)
        np.testing.assert_allclose(points[:, 1].max() - points[:, 1].min(), 2)

    def test_edge_touching_mask_closes(self):
        mask = decode_voxel_mask(encoded([[0, 4]]), base_width=2, base_height=2, xy_scale=1)
        image = ImageVolume(np.zeros((1, 2, 2), dtype='u2'), np.eye(4), np.eye(4), (), 'CT', 1, -1024, '')
        self.assertEqual(len(mask_to_contours(mask, image)), 1)

    def test_registered_planes_use_individual_affines(self):
        mask = decode_voxel_mask(encoded([[0, 1, 3], [0, 1, 3]]), base_width=2, base_height=2, xy_scale=1)
        first, second = np.eye(4), np.eye(4)
        second[:3, 3] = [4, -3, 10]
        image = ImageVolume(np.zeros((2, 2, 2), dtype='u2'), np.eye(4), np.eye(4), (),
                            'CT', 1, -1024, '', (first, second), False)
        contours = mask_to_contours(mask, image)
        np.testing.assert_allclose(contours[1][1].mean(axis=0) - contours[0][1].mean(axis=0), [4, -3, 10])


if __name__ == '__main__':
    unittest.main()
