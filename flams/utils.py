#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: kasgel
"""

from pathlib import Path
import appdirs

""" setup
This script deals with getting a platform-specific directory to store app data.
E.g. on     Linux: ~/.local/share/<AppName>
            Mac: '/Users/trentm/Library/Application Support/SuperApp'
            Windows: 'C:\\Users\\trentm\\AppData\\Local\\Acme\\SuperApp'
"""

def get_data_dir(app_name="flams"):
    """
    This function gets a platform-specific directory to store app data.

    Parameters
    ----------
    app_name: str
    Name of application, i.e., flams

    """
    # Ensure data dir exists and return.
    data_dir = appdirs.user_data_dir(app_name)
    Path(data_dir).mkdir(parents=True, exist_ok=True)
    return data_dir
