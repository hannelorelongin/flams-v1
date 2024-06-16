import unittest
import tempfile
import os
from types import SimpleNamespace
from pathlib import Path

from src.flams.input import check_data_dir

class InputTest(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_check_data_dir(self):
        # Nonexistent directory, creates it:
        path = Path(f"{self.temp_dir.name}/nonexistent")
        self.assertFalse(path.exists())
        args = SimpleNamespace(data_dir=path)
        check_data_dir(args)
        self.assertTrue(path.exists())

        # Existing directory, no error:
        path = Path(f"{self.temp_dir.name}/existing_dir")
        os.mkdir(path)
        args = SimpleNamespace(data_dir=path)
        check_data_dir(args)

        # Existing file, error:
        path = Path(f"{self.temp_dir.name}/existing_file")
        path.touch()

        args = SimpleNamespace(data_dir = path)
        self.assertRaises(SystemExit, check_data_dir, args)

if __name__ == '__main__':
    unittest.main()
