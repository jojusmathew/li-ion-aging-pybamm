"""Download the per-cell summary tables of the LG M50T ageing study without the 6 GB zip.

Data: Kirkaldy et al., J. Power Sources 603 (2024) 234185, Zenodo 10.5281/zenodo.10637534 (CC-BY-4.0).
Only experiment 2,2 is used: cycling between 70 and 85 % SoC at 10, 25 and 40 degC, two cells each.
A zip keeps its table of contents at the end, so HTTP range requests are enough to pull single files.
"""

import io
import urllib.request
import zipfile
from pathlib import Path

URL = "https://zenodo.org/records/10637534/files/Expt%202%2C2%20-%20C-based%20Degradation%202.zip?download=1"
FOLDER = "Expt 2,2 - C-based Degradation 2/Summary Data/Performance Summary/"


class RemoteFile(io.RawIOBase):
    """Read-only, seekable view of a remote file, fetched in pieces with HTTP range requests."""

    def __init__(self, url):
        self.url, self.pos = url, 0
        head = urllib.request.urlopen(urllib.request.Request(url, method="HEAD"))
        self.size = int(head.headers["Content-Length"])

    def readable(self):
        return True

    def seekable(self):
        return True

    def tell(self):
        return self.pos

    def seek(self, offset, whence=io.SEEK_SET):
        self.pos = (offset, self.pos + offset, self.size + offset)[whence]
        return self.pos

    def readinto(self, buffer):
        if self.pos >= self.size:
            return 0
        end = min(self.pos + len(buffer), self.size) - 1
        request = urllib.request.Request(self.url, headers={"Range": f"bytes={self.pos}-{end}"})
        data = urllib.request.urlopen(request).read()
        buffer[: len(data)] = data
        self.pos += len(data)
        return len(data)


if __name__ == "__main__":
    out = Path(__file__).parent / "data"
    out.mkdir(exist_ok=True)
    with zipfile.ZipFile(io.BufferedReader(RemoteFile(URL), buffer_size=1 << 16)) as z:
        for name in z.namelist():
            if name.startswith(FOLDER) and name.endswith(".csv"):
                (out / Path(name).name).write_bytes(z.read(name))
                print("saved", Path(name).name)
