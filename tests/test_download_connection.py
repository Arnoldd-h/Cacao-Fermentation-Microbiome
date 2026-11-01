"""A local stalled HTTP stream must time out, preserve a partial and resume."""

import hashlib
import importlib.util
import sys
import tempfile
import threading
import unittest
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
from cacao_inventory.download import download_file

spec = importlib.util.spec_from_file_location("download_read", ROOT / "scripts/qc/download_pilot_read.py")
DOWNLOAD = importlib.util.module_from_spec(spec)
spec.loader.exec_module(DOWNLOAD)


class DownloadConnectionTests(unittest.TestCase):
    def test_stalled_transfer_keeps_partial_and_retry_resumes_by_range(self):
        release = threading.Event()
        requests = []

        class Handler(BaseHTTPRequestHandler):
            def log_message(self, *args):
                pass

            def do_GET(self):
                offset = self.headers.get("Range", "")
                requests.append(offset)
                self.send_response(206 if offset else 200)
                self.send_header("Content-Length", "3" if offset else "4")
                self.end_headers()
                if offset:
                    self.wfile.write(b"BCD")
                else:
                    self.wfile.write(b"A")
                    self.wfile.flush()
                    release.wait(5)

        server = ThreadingHTTPServer(("127.0.0.1", 0), Handler)
        worker = threading.Thread(target=server.serve_forever, daemon=True)
        worker.start()
        try:
            with tempfile.TemporaryDirectory() as temporary:
                target = Path(temporary) / "synthetic.fastq.gz"
                checksum = hashlib.md5(b"ABCD").hexdigest()
                opener = DOWNLOAD.download_opener({"inventory": {"http_timeout_seconds": 0.2}})
                url = "http://127.0.0.1:" + str(server.server_port) + "/fixture"
                with self.assertRaises(TimeoutError):
                    download_file(url, target, 4, checksum, 1, opener=opener)
                self.assertFalse(target.exists())
                self.assertEqual(target.with_name(target.name + ".part").read_bytes(), b"A")
                release.set()
                self.assertEqual(download_file(url, target, 4, checksum, 1, opener=opener), "resumed_valid")
                self.assertEqual(target.read_bytes(), b"ABCD")
                self.assertEqual(requests, ["", "bytes=1-"])
        finally:
            release.set()
            server.shutdown()
            server.server_close()
            worker.join(timeout=5)

    def test_invalid_socket_timeout_fails_before_download(self):
        for value in (0, -1, float("nan"), float("inf")):
            with self.subTest(value=value), self.assertRaises(ValueError):
                DOWNLOAD.download_opener({"inventory": {"http_timeout_seconds": value}})


if __name__ == "__main__":
    unittest.main()
