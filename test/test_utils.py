from pathlib import Path
import unittest
from flams.utils import get_data_dir
import appdirs

APPDIR_FOR_TESTS = "flams-unittest"


class UtilsTestCase(unittest.TestCase):
    def tearDown(self) -> None:
        data_dir = appdirs.user_data_dir(APPDIR_FOR_TESTS)
        Path(data_dir).rmdir()

    def test_data_dir(self):
        data_dir = appdirs.user_data_dir(APPDIR_FOR_TESTS)

        # Assert that the data dir does NOT exist
        self.assertFalse(Path(data_dir).exists())

        # Call get_data_dir which will also create the folder if it does not exist
        get_data_dir(app_name=APPDIR_FOR_TESTS)

        # Assert that the data dir now DOES exist.
        self.assertTrue(Path(data_dir).exists())
