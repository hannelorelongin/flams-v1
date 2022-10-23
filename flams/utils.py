from pathlib import Path
import appdirs


# This will get a platform-specific dir to store app data.
# E.g. on Linux: ~/.local/share/<AppName>
# Mac: '/Users/trentm/Library/Application Support/SuperApp'
# Windows: 'C:\\Users\\trentm\\AppData\\Local\\Acme\\SuperApp'
def get_data_dir():
    # Ensure data dir exists and return.
    data_dir = appdirs.user_data_dir("flams")
    Path(data_dir).mkdir(parents=True, exist_ok=True)
    return data_dir
