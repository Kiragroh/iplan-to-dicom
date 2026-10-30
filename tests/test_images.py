import io
import unittest
import xml.etree.ElementTree as ET
from pathlib import Path
import tempfile

import numpy as np
from PIL import Image
from pydicom.dataset import Dataset

from iplan_dicom.images import ImageVolume, checked_center_matrix, decode_png, registration_matrix, regularize_image_volume
from iplan_dicom.xbrain import FormatError


class ImageTests(unittest.TestCase):
    def test_regularization_identity_preserves_all_samples(self):
        pixels = (np.arange(3 * 4 * 5).reshape(3, 4, 5) + 100).astype('u2')
        image = ImageVolume(pixels, np.eye(4), np.eye(4), tuple(Dataset() for _ in range(3)), 'CT', 1, -1024, 'synthetic')
        result = regularize_image_volume(image, np.eye(4))
        np.testing.assert_array_equal(result.pixels, pixels)
        np.testing.assert_allclose(result.voxel_to_patient, np.eye(4))
        self.assertTrue(result.uniform_grid)

    def test_regularization_interpolates_physical_plane_positions(self):
        yy, xx = np.indices((5, 5))
        planes, matrices = [], []
        for dx, z in [(0, 0), (.2, 1.2), (-.2, 2.4)]:
            matrix = np.eye(4)
            matrix[:3, 3] = [dx, 0, z]
            matrices.append(matrix)
            planes.append(np.rint(100 + 10 * (xx + dx) + 20 * yy + 30 * z).astype('u2'))
        image = ImageVolume(np.stack(planes), np.eye(4), np.eye(4), tuple(Dataset() for _ in range(3)),
                            'CT', 1, -1024, 'synthetic', tuple(matrices), False)
        result = regularize_image_volume(image, np.eye(4))
        for z in range(3):
            for y in range(1, 4):
                for x in range(1, 4):
                    self.assertEqual(result.pixels[z, y, x + 1], 100 + 10 * x + 20 * y + 30 * z)
        self.assertTrue(np.all(result.pixels[-1] == 0))  # outside source coverage = air
        np.testing.assert_array_equal(image.pixels, np.stack(planes))

    def test_regularization_rejects_crossing_source_planes(self):
        first, second = np.eye(4), np.eye(4)
        angle = .1
        second[:3, :3] = [[np.cos(angle), 0, np.sin(angle)], [0, 1, 0], [-np.sin(angle), 0, np.cos(angle)]]
        second[:3, 3] = [0, 0, .1]
        image = ImageVolume(np.full((2, 5, 5), 1024, dtype='u2'), np.eye(4), np.eye(4), (Dataset(), Dataset()),
                            'CT', 1, -1024, 'synthetic', (first, second), False)
        with self.assertRaises(FormatError):
            regularize_image_volume(image, np.eye(4))

    def fusion_document(self, edges):
        root = ET.Element('File')
        version = ET.SubElement(ET.SubElement(root, 'Content'), 'ContentVersion', value='1')
        for source, target, translation in edges:
            fusion = ET.SubElement(version, 'Fusion')
            ET.SubElement(fusion, 'Reference', contentType='SliceSet', contentUid=source)
            ET.SubElement(fusion, 'Reference', contentType='SliceSet', contentUid=target)
            matrix = np.eye(4)
            matrix[:3, 3] = translation
            value = ET.SubElement(ET.SubElement(ET.SubElement(ET.SubElement(fusion, 'Transformation'), 'Matrix_44'), '_DOUBLEVEC_'), '_ASCII_')
            value.text = ','.join(map(str, [16] + matrix.ravel().tolist()))
        return ET.tostring(root)

    def test_fusion_forward_chain_and_direct_precedence(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'fusion.xbrain'
            edges = [('a', 'b', [10, 0, 0]), ('b', 'c', [0, 4, 0])]
            path.write_bytes(self.fusion_document(edges))
            np.testing.assert_allclose(registration_matrix('a', 'c', Path(folder))[:3, 3], [10, 4, 0])
            np.testing.assert_allclose(registration_matrix('c', 'a', Path(folder))[:3, 3], [-10, -4, 0])
            path.write_bytes(self.fusion_document(edges + [('a', 'c', [9, 4, 0])]))
            np.testing.assert_allclose(registration_matrix('a', 'c', Path(folder))[:3, 3], [9, 4, 0])

    def test_fusion_missing_and_conflicting_paths_fail(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'fusion.xbrain'
            path.write_bytes(self.fusion_document([('a', 'b', [1, 0, 0]), ('a', 'b', [2, 0, 0])]))
            with self.assertRaises(FormatError):
                registration_matrix('a', 'b', Path(folder))
            with self.assertRaises(FormatError):
                registration_matrix('a', 'missing', Path(folder))

    def test_brainlab_png_byte_order_and_stored_range(self):
        pixels = np.array([[0, 1, 256], [1024, 4095, 27]], dtype=np.uint16)
        stream = io.BytesIO()
        Image.fromarray(pixels.byteswap()).save(stream, format="PNG")
        actual = decode_png(stream.getvalue(), bits=16, bits_stored=12)
        np.testing.assert_array_equal(actual, pixels)
        with self.assertRaises(FormatError):
            decode_png(stream.getvalue(), bits=16, bits_stored=10)

    def test_eight_bit_png_is_not_swapped(self):
        pixels = np.array([[1, 255], [0, 13]], dtype=np.uint8)
        stream = io.BytesIO()
        Image.fromarray(pixels).save(stream, format="PNG")
        np.testing.assert_array_equal(decode_png(stream.getvalue(), bits=8), pixels)

    def test_voxel_corners_and_left_right_reflection(self):
        # Four x pixels: original DICOM centers 10,12,14,16 mm.
        native = np.array([[2, 0, 0, -17], [0, -3, 0, -18.5], [0, 0, -1, 30], [0, 0, 0, 1]])
        values = [18] + native.ravel().tolist() + [2, -1]
        node = ET.fromstring('<Image width="4" height="3" headFeetOrientation="HeadFirst" patientOrientation="Supine"><SubSystem><Point_3D_Subsystem><_DOUBLEVEC_><_ASCII_>' + ','.join(map(str, values)) + '</_ASCII_></_DOUBLEVEC_></Point_3D_Subsystem></SubSystem></Image>')
        header = Dataset()
        header.PatientPosition = 'HFS'
        header.Rows, header.Columns = 3, 4
        header.ImagePositionPatient = [10, 20, 30]
        header.ImageOrientationPatient = [1, 0, 0, 0, 1, 0]
        header.PixelSpacing = [3, 2]
        actual = checked_center_matrix(node, header)
        np.testing.assert_allclose(actual @ [0, 0, 0, 1], [16, 20, 30, 1])
        np.testing.assert_allclose(actual @ [3, 2, 0, 1], [10, 26, 30, 1])
        header.ImagePositionPatient = [11, 20, 30]
        with self.assertRaises(FormatError):
            checked_center_matrix(node, header)


if __name__ == '__main__':
    unittest.main()
